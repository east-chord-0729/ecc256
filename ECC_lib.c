#include "ECC_lib.h"

// E/K: y^2 = x^3 - 3x + b
// p = 2^256 − 2^224 + 2^192 + 2^96 − 1
// jacobian projective: c = 2, d = 3

// deep copy of ec_point_af and ec_point_pj
void set_ec_point_af(EC_POINT_AF *dest, const EC_POINT_AF* src)
{
    set_bn(&dest->x, &src->x);
    set_bn(&dest->y, &src->y);
    dest->is_infty = src->is_infty;
}
void set_ec_point_pj(EC_POINT_PJ *dest, const EC_POINT_PJ* src)
{
    set_bn(&dest->x, &src->x);
    set_bn(&dest->y, &src->y);
    set_bn(&dest->z, &src->z);
    dest->is_infty = src->is_infty;
}

// jaco to affn: (X:Y:Z) --> (X/Z^2, Y/Z^3)
static void jc2af(EC_POINT_AF* point_af, const EC_POINT_PJ* point_pj)
{
    BN inv_z2;      // inverse of Z^2
    BN inv_z3;      // inverse of Z^3

    // check infty
    if(point_pj->is_infty) {
        point_af->is_infty = 1;
        return;
    }
    point_af->is_infty = 0;

    // inverse of Z^2 and Z^3
    sqrp(&inv_z2, &point_pj->z);
    mulp(&inv_z3, &inv_z2, &point_pj->z);
    inv(&inv_z2, &inv_z2);
    inv(&inv_z3, &inv_z3);

    // jaco to affn: x = X/Z^2
    mulp(&point_af->x, &point_pj->x, &inv_z2);
    
    // jaco to affn: y = Y/Z^3
    mulp(&point_af->y, &point_pj->y, &inv_z3);
}

// affine to jacobian : (x, y) --> (x:y:1)
static void af2jc(EC_POINT_PJ* point_pj, const EC_POINT_AF* point_af)
{
    // check infty
    if(point_af->is_infty) {
        point_pj->is_infty = 1;
        return;
    }
    point_pj->is_infty = 0;

    // affn to jaco : X = x
    point_pj->x = point_af->x;

    // affn to jaco : Y = y
    point_pj->y = point_af->y;

    // affn to jaco : Z = 1
    point_pj->z.v[0] = 1;
    for(int i = 1; i < NUMWORD; i++) {
        point_pj->z.v[i] = 0;
    }
}

// addtion over affine, 1-inv + 2-mul + 1-sqr + 6-add
static void ecadd_af(EC_POINT_AF* point_r, const EC_POINT_AF* point_p, const EC_POINT_AF* point_q) 
{
    BN lambda, tmp;
    BN rx, ry;

    // check infty
    if (point_p->is_infty) {
        set_ec_point_af(point_r, point_q);
        return;
    } else if (point_q->is_infty) {
        set_ec_point_af(point_r, point_p);
        return;
    }

    // check inverse available
    if (!ucmp(&point_p->x, &point_q->x)) return;

    // lambda = (Qy - Py) / (Qx - Px)
    subp(&lambda, &point_q->y, &point_p->y);
    subp(&tmp, &point_q->x, &point_p->x);
    inv(&tmp, &tmp);
    mulp(&lambda, &lambda, &tmp);

    // Rx = lambda^2 - Px - Qx
    sqrp(&rx, &lambda);
    subp(&rx, &rx, &point_p->x);
    subp(&rx, &rx, &point_q->x);
 
    // Ry = L * (Rx - Px) - Py
    subp(&ry, &point_p->x, &rx);
    mulp(&ry, &ry, &lambda);
    subp(&ry, &ry, &point_p->y);

    // return, *note: ecadd_af(&R, &R, &Q)
    set_bn(&point_r->x, &rx);
    set_bn(&point_r->y, &ry);
    point_r->is_infty = 0;
}

// doubling over affine, 1-inv + 2-mul + 2-sqr + 7-add
static void ecdbl_af(EC_POINT_AF* point_r, const EC_POINT_AF* point_p) 
{
    BN lambda, tmp;
    BN rx, ry;

    // check infty
    if (point_p->is_infty) {
        point_r->is_infty = 1;
        return;
    }

    // check inverse available
    if (!ucmp(&point_p->y, &zero)) return;

    // lambda = (3 * Px^2 + a) / 2Py
    sqrp(&tmp, &point_p->x);
    addp(&lambda, &tmp, &tmp);
    addp(&lambda, &lambda, &tmp);
    addp(&lambda, &lambda, &coef_a);
    addp(&tmp, &point_p->y, &point_p->y);
    inv(&tmp, &tmp);
    mulp(&lambda, &lambda, &tmp);

    // Rx = L^2 - 2Px
    sqrp(&rx, &lambda);
    subp(&rx, &rx, &point_p->x);
    subp(&rx, &rx, &point_p->x);

    // Ry = L * (Px - Rx) - Py
    subp(&ry, &point_p->x, &rx);
    mulp(&ry, &ry, &lambda);
    subp(&ry, &ry, &point_p->y);

    // return, *note: ecdbl_af(&R, &R, &Q)
    set_bn(&point_r->x, &rx);
    set_bn(&point_r->y, &ry);
    point_r->is_infty = 0;
}

/* 
    아핀에서는 역원을 사용하기 때문에 느리다. 하지만 사영좌표계를 이용하면 역원 연산을 쓰지 않는다.
    아핀좌표계를 사영좌표계로 바꿔준 후 점연산을 진행한다. 연산을 마치면 다시 아핀좌표계로 바꾸어준다.
    ( af2jc --> scalar multiplication --> jc2af )
*/

// doubling over jacobian, jacobian = 2*jacobian, 4-mul + 4sqr + 9-add
static void ecdbl_jc(EC_POINT_PJ* point_r, const EC_POINT_PJ* point_p) 
{
    BN t1, t2, t3;
    BN rx, ry, rz;

    // check infty
    if (point_p->is_infty) {
        point_r->is_infty = 1;
        return;
    }

    // Guide to ECC, p.91
    sqrp(&t1, &point_p->z);
    subp(&t2, &point_p->x, &t1);
    addp(&t1, &point_p->x, &t1);
    mulp(&t2, &t2, &t1);
    addp(&t3, &t2, &t2);
    addp(&t2, &t2, &t3);
    addp(&ry, &point_p->y, &point_p->y);
    mulp(&rz, &ry, &point_p->z);
    sqrp(&ry, &ry);
    mulp(&t3, &ry, &point_p->x);
    sqrp(&ry, &ry);
    
     // shift, (ry >> 1) < p
    if(ry.v[0] & 1) { 
        // if ry is odd
        uint32_t carry = uadd(&ry, &ry, &P);
        rshift1(&ry, &ry);
        ry.v[NUMWORD-1] += (carry << (WORDBITS-1));
    } else {
        // if ry is even
        rshift1(&ry, &ry);
    }

    sqrp(&rx, &t2);
    addp(&t1, &t3, &t3);
    subp(&rx, &rx, &t1);
    subp(&t1, &t3, &rx);
    mulp(&t1, &t1, &t2);
    subp(&ry, &t1, &ry);

    // return, *note: ecdbl_jc(&R, &R)
    set_bn(&point_r->x, &rx);
    set_bn(&point_r->y, &ry);
    set_bn(&point_r->z, &rz);
    point_r->is_infty = 0;
}

// addtion over jacobian, jacobian = jacobian + affine, 8-mul + 3-sqr + 7-add
static void ecadd_jc(EC_POINT_PJ* point_r, const EC_POINT_PJ* point_p, const EC_POINT_AF* point_q) 
{
    BN t1, t2, t3, t4;
    BN rx, ry, rz;

    // check infty
    if (point_p->is_infty) {
        af2jc(point_r, point_q);
        return;
    } else if (point_q->is_infty) {
        set_ec_point_pj(point_r, point_p);
        return;
    }

    // Guide to ECC, p.91
    sqrp(&t1, &point_p->z);
    mulp(&t2, &t1, &point_p->z);
    mulp(&t1, &t1, &point_q->x);
    mulp(&t2, &t2, &point_q->y);
    subp(&t1, &t1, &point_p->x);
    subp(&t2, &t2, &point_p->y);

    // t1 = Pz^2 * (Qx - Px)
    if (!ucmp(&t1, &zero)) 
    {
        // t2 = Pz^3 * (Qy - Py)
        if (!ucmp(&t2, &zero)) 
        {
            // t1, t2 = 0 --> Qx = Px and Qy = Py
            af2jc(point_r, point_q);
            ecdbl_jc(point_r, point_r);
            return; 
        } 
        else {
            // t1 = 0 --> Rz = t1 * Pz = 0
            point_r->is_infty = 1;
            return;
        }
    }

    // Guide to ECC, p.92
    mulp(&rz, &point_p->z, &t1);
    sqrp(&t3, &t1);
    mulp(&t4, &t3, &t1);
    mulp(&t3, &t3, &point_p->x);
    addp(&t1, &t3, &t3);
    sqrp(&rx, &t2);
    subp(&rx, &rx, &t1);
    subp(&rx, &rx, &t4);
    subp(&t3, &t3, &rx);
    mulp(&t3, &t3, &t2);
    mulp(&t4, &t4, &point_p->y);
    subp(&ry, &t3, &t4);

    // return, *note: ecadd_jc(&R, &R, &Q)
    set_bn(&point_r->x, &rx);
    set_bn(&point_r->y, &ry);
    set_bn(&point_r->z, &rz);
    point_r->is_infty = 0;
}

/*  
    general scalar multiplication: O(2^n)
    LtoR, RtoL and so on: O(n)  
    LtoR, RtoL: non-constant, Montgomery ladder: constant
*/

// scalar multiplication of ec, left to right
void ecsm_ltr(EC_POINT_AF* point_r, const EC_POINT_AF* point_G, const BN* scalar)
{
    EC_POINT_PJ ret_pj = {0};

    // init
    ret_pj.is_infty = 1;

    // left to right algorithm
    for(int i = NUMWORD - 1; i >= 0; i--) 
    {
        for (int j = WORDBITS - 1; j >= 0; j--) 
        {
            // jacobian = 2 * jacobian
            ecdbl_jc(&ret_pj, &ret_pj);

            if ((scalar->v[i] >> j) & 1) {
                // jacobian = jacobian + affine
                ecadd_jc(&ret_pj, &ret_pj, point_G);
            }
        }
    }

    // jacobian to affine
    jc2af(point_r, &ret_pj);
}

// scalar multiplication of ec, left to right
void ecsm_rtl(EC_POINT_AF* point_r, const EC_POINT_AF* point_G, const BN* scalar)
{
    EC_POINT_AF g_af = {0};
    EC_POINT_PJ ret_pj = {0};

    // init
    set_ec_point_af(&g_af, point_G);
    ret_pj.is_infty = 1;

    // right to left algorithm
    for(int i = 0; i < NUMWORD; i++) 
    {
        for (int j = 0; j < WORDBITS; j++) 
        {
            if ((scalar->v[i] >> j) & 1) {
                // jacobian = jacobian + affine
                ecadd_jc(&ret_pj, &ret_pj, &g_af);
            }

            // affine = affine + affine //! use jacobian only coordinates
            ecdbl_af(&g_af, &g_af);
        }                  
    }

    // proj --> affn
    jc2af(point_r, &ret_pj);
}

/*  
    left to right algorithm pre-computed version
    scalar를 상위 n비트씩 한번에 읽어서 곱셈을 한번에 처리할 수 있다. (여기서는 8비트씩 읽음)

    right to left algorithm pre-computed version
    더블링하는 point_G가 고정 값이므로 사전 계산하여 연산 가능.
*/

// left to right algorithm pre-computed version
void ecsm_ltr_precomp(EC_POINT_AF* point_r, const EC_POINT_AF* point_G, const BN* scalar)
{
    uint8_t offset = 0;
    EC_POINT_PJ ret_pj = {0};

    // init
    ret_pj.is_infty = 1;

    // left to right algorithm
    for(int i = NUMWORD - 1; i >= 0; i--) 
    {
        for (int k = 3; k >= 0; k--) 
        {
            for (int j = 0; j < 8; j++){
                // doubling
                ecdbl_jc(&ret_pj, &ret_pj);
            }

            // scanning 8-bits
            offset = (scalar->v[i] >> (k * 8)) & 0xFF;

            // addition
            ecadd_jc(&ret_pj, &ret_pj, &fix_g_ltr[offset]);
        }
    }

    // proj --> affn
    jc2af(point_r, &ret_pj);
}

// right to left algorithm pre-computed version
void ecsm_rtl_precomp(EC_POINT_AF* point_r, const EC_POINT_AF* point_G, const BN* scalar)
{
    EC_POINT_PJ ret_pj = {0};

    // init
    ret_pj.is_infty = 1;

    // right to left algorithm
    for(int i = 0; i < NUMWORD; i++) 
    {
        for (int j = 0; j < WORDBITS; j++) 
        {
            if ((scalar->v[i] >> j) & 1) {
                // addition
                ecadd_jc(&ret_pj, &ret_pj, &fixG_RtoL[i*32+j]);
            }
        }
    }

    // proj --> affn
    jc2af(point_r, &ret_pj);
}





/*  LtoR, RtoL은 non-constant --> 몽고메리래더 알고리즘을 사용. 
    몽고메리 래더는 병렬연산이 가능하다는 장점이 있다. */
// void ecsm_mont(EC_POINT_AF *point_r, const EC_POINT_AF *point_g, const BN* scalar) 
// {
//     EC_POINT_PJ left = {0};
//     EC_POINT_PJ right = {0};

//     // init : left = infty, right = G
//     af2jc(&right, point_g);

//     // run montgomery ladder
//     for(int i = NUMWORD - 1; i >= 0; i--) 
//     {
//         for (int j = WORDBITS - 1; j >= 0; j--) 
//         {
//             if ((scalar->v[i] >> j) & 1) 
//             {
//                 ecadd_af(&left, &left, &right); //? j = j + a --> a가 아핀이 아닌데요?
//                 ecdbl_jc(&right, &right);
//             } 
//             else 
//             {
//                 ecadd_jc(&right, &right, &left);
//                 ecdbl_jc(&left, &left);
//             }
//         }
//     }

//     /*
//         mix 좌표계가 안되는데 어떻게 할까?
//         1. 될수있는 기가막힌 아이디어가 있을 수 있다.
//         2. 그냥 mix 없이 자코비안 좌표 혹은 아핀 좌표만을 사용한다. 
//     */
// }

int main(void)
{
    const BN k = {
        0x3933224B, 0x18671BCA, 0x5E4D9E0A, 0xBA08EE99, 
        0xB568A7A2, 0xB6D14865, 0x71AFC9F6, 0xDDB7F114};

    const EC_POINT_AF G = {
        {0xD898C296, 0xF4A13945, 0x2DEB33A0, 0x77037D81,
         0x63A440F2, 0xF8BCE6E5, 0xE12C4247, 0x6B17D1F2},
        {0x37BF51F5, 0xCBB64068, 0x6B315ECE, 0x2BCE3357,
         0x7C0F9E16, 0x8EE7EB4A, 0xFE1A7F9B, 0x4FE342E2}, 0};

    EC_POINT_AF R = {0};
    
    ecsm_ltr(&R, &G, &k);

    for(int i = 0; i < 8; i++) printf("[%d]%08x ", i, R.x.v[i]); puts("");
    for(int i = 0; i < 8; i++) printf("[%d]%08x ", i, R.y.v[i]); puts("");

    return 0;
}

/*

#1

0. R을 무한원점으로 설정한다.
1. R을 더블링한다.
2. k를 앞 비트부터 읽는다.
3. 1이면 P와 덧셈을 하고, 0이면 하지 않는다.

k = 10110_(2진수) = 22_(10진수) --> k4 = 1, k3 = 0, k2 = 1, k1 = 1, k0 = 0

R = infty

k4 = 1
    R = 2R = infty
    R = R+P = P
k3 = 0
    R = 2R = 2P
k2 = 1
    R = 2R = 4P
    R = R+P = 5P
k1 = 1
    R = 2R = 10P
    R = R+P = 11P
k0 = 0
    R = 2R = 22P --> k = 22 로 계산 정답!
*/





/*

    EC_POINT precomputed_point;
    set_ec_point(&precomputed_point, &G);

    for(int i = 0; i < 256; i++) {
        printf("{{");
        for(int j = 0; j < 7; j++) printf("0x%08x, ", precomputed_point.x.v[j]);
        printf("0x%08x}, ", precomputed_point.x.v[7]);

        printf("\n");

        printf("{");
        for(int j = 0; j < 7; j++) printf("0x%08x, ", precomputed_point.y.v[j]);
        printf("0x%08x}, ", precomputed_point.y.v[7]);

        printf("\n");

        printf("{1,0,0,0,0,0,0,0}}, ");

        printf("\n");

        ecdbl_af(&precomputed_point, &precomputed_point);
    }

*/