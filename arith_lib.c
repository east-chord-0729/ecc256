#include "arith_lib.h"

const BN nistp256[5] = {
    {0xffffffff, 0xffffffff, 0xffffffff, 0x00000000, 
     0x00000000, 0x00000000, 0x00000001, 0xffffffff},
    {0xfffffffe, 0xffffffff, 0xffffffff, 0x00000001, 
     0x00000000, 0x00000000, 0x00000002, 0xfffffffe},
    {0xfffffffd, 0xffffffff, 0xffffffff, 0x00000002,
     0x00000000, 0x00000000, 0x00000003, 0xfffffffd},
    {0xfffffffc, 0xffffffff, 0xffffffff, 0x00000003,
     0x00000000, 0x00000000, 0x00000004, 0xfffffffc},
    {0xfffffffb, 0xffffffff, 0xffffffff, 0x00000004,
     0x00000000d, 0x00000000, 0x00000005, 0xfffffffb}};

const BN P_prime = {
    0x00000001, 0x00000000, 0x00000000, 0x00000001, 
    0x00000000, 0x00000000, 0x00000002, 0xFFFFFFFF};

const BN RmodP = {
    0x00000003, 0x00000000, 0xFFFFFFFF ,0xFFFFFFFB, 
    0xFFFFFFFE, 0xFFFFFFFF, 0xFFFFFFFD, 0x00000004};

const BN RRmodP = {
    0x00000003, 0x00000000, 0xFFFFFFFF ,0xFFFFFFFB, 
    0xFFFFFFFE, 0xFFFFFFFF, 0xFFFFFFFD, 0x00000004};

/* BN deep copy */
//todo void set_bn(BN* dest, const int32_t* src)
//todo void _bn(BN* dest, const int32_t* src)
void set_bn(BN* dest, const BN* src)
{
    for (int i = 0; i < NUMWORD; i++){
        dest->v[i] = src->v[i];
    }
    dest->s = src->s;
}

/* compare function: memcmp는 하위 워드부터 비교하므로 따로 구현함 */
int32_t ucmp(const BN* opa, const BN* opb) {
    for (int i = NUMWORD - 1; i >= 0; i--) {
        if (opa->v[i] > opb->v[i]) return 1;        // opa > opb
        if (opa->v[i] < opb->v[i]) return -1;       // opa < opb
    }
    return 0;   // opa = opb
}

/* ret = opa >> 1 */
uint32_t rshift1(BN* ret, const BN* opa)
{
    BN r;
    uint32_t carry = 0;

    for (int i = NUMWORD - 1; i >= 0; i--) {
        r.v[i] = (opa->v[i] >> 1) ^ carry;       // 1 right shift
        carry = (opa->v[i] & 1) << 31;           // carry
    }

    set_bn(ret, &r);

    return carry;
} 

/* ret = opa + opb */
uint32_t uadd(BN* ret, const BN* opa, const BN* opb)
{
    BN r = {0};
    uint32_t carry = 0;

    for(size_t i = 0; i < NUMWORD; i++) {
        r.v[i] = opa->v[i] + opb->v[i] + carry;                  // addition
        carry = (opa->v[i] + opb->v[i] + carry < carry) +       // carry condition 1
              (opa->v[i] + opb->v[i] < opa->v[i]);              // carry condition 2
    }

    set_bn(ret, &r);
 
    return carry;
}

/* ret = opa - opb , return carry */
uint32_t usub(BN* ret, const BN* opa, const BN* opb)
{
    BN r;
    uint32_t carry = 0;

    for (int i = 0; i < NUMWORD; i++) {
        r.v[i] = opa->v[i] - opb->v[i] - carry;                 // substraction
        carry = (opa->v[i] < opb->v[i]) ||                      // carry condition 1
              ((opa->v[i] == opb->v[i]) && (carry == 1));       // carry condition 2
    }
    
    set_bn(ret, &r);

    return carry;
}

/* ret = opa + opb mod p */
void addp(BN* ret, const BN* opa, const BN* opb) 
{
    BN r;
    uint32_t carry = 0;

    carry = uadd(&r, opa, opb);

    if (carry == 1 || ucmp(&r, &P) >= 0) {
        usub(&r, &r, &P);
    }

    set_bn(ret, &r);
}

/* ret = opa - opb mod p */
void subp(BN *ret, const BN* opa, const BN* opb) 
{
    BN r;
    uint32_t carry = 0;

    carry = usub(&r, opa, opb);

    if (carry == 1) {
        uadd(&r, &r, &P);
    }

    set_bn(ret, &r);
}

/*  R = A^2 product scanning
    mul_ps랑 같은 방식으로 동작하되, 제곱의 경우 uv값이 한 번 중복되어 더해지는 경우가 있음.
    ex. uv = a2*a0 + a1*a1 + a0*a2 의 경우 a0*a2 = a2*a0이므로 한 번만 곱셈 후 시프트. */
void usqr_ps(BN2* ret, const BN* opa)
{
    uint32_t r2 = 0, r1 = 0, r0 = 0, c1 = 0, c2 = 0;
    uint64_t uv = 0;
    uint32_t u = 0, v = 0;
    
    for (int k = 0; k < 2*NUMWORD-1; k++) 
    {
        for(int i = 0, j = k; i <= k; i++, j--) 
        {
            if (i >= NUMWORD || j >= NUMWORD || j < i) continue;

            // uv 계산. 중복되는 부분은 시프트로.
            uv = (uint64_t)opa->v[i] * opa->v[j];
            if (i < j) {
                c2 = (uint32_t)(uv >> 63);
                uv = uv << 1;
                r2 += c2;
            }
            v = (uint32_t)uv;
            u = (uint32_t)(uv >> 32);

            // uv를 계속 r2||r1||r0에 담으면서 더함.
            r0 += v;
            c1 = (r0 < v);
            r1 += u;
            c2 = (r1 < u);
            r1 += c1;
            c2 += (r1 < c1);
            r2 += c2;
        }

        // 한칸 올라감.
        ret->v[k] = r0;
        r0 = r1;
        r1 = r2;
        r2 = 0;
    }
    ret->v[2*NUMWORD-1] = r0;
}

/*  R = A * B school book version: 우리가 일반적으로 하는 곱셈과 똑같은 방식으로 진행됨.
    32비트 x 32비트는 64비트 자료형에 담고 상위32비트와 하위32비트로 나눔.
    자세한 내용은 하단의 #1을 참고. */
void umul_os(BN2* ret, const BN* opa, const BN* opb)
{                               
    uint64_t uv = 0;
    uint32_t u = 0;

    // uv = ret->v[i+j] + (uint64_t)opa->v[i] * opb->v[j] + u; 부분이 동작하기 위해 초기화.
    memset(ret, 0, sizeof(BN2));
    
    // A * B -- 하단 #1 참고
    for(int i = 0; i < NUMWORD; i++) {
        u = 0;
        for(int j = 0; j < NUMWORD; j++) {
            uv = ret->v[i+j] + (uint64_t)opa->v[i] * opb->v[j] + u;
            ret->v[i+j] = (uint32_t)uv;
            u = (uint32_t)(uv >> 32);
        }
        ret->v[i+NUMWORD] = u;
    }
}

/*  R = A * B product scanning: 우리가 일반적으로 하는 곱셈과 똑같은 방식으로 진행됨.  
    자세한 내용은 하단의 #2을 참고. */
void umul_ps(BN2* ret, const BN* opa, const BN* opb)
{
    uint32_t r2 = 0, r1 = 0, r0 = 0, c1 = 0, c2 = 0;
    uint64_t uv = 0;
    uint32_t u = 0, v = 0;

    for (int k = 0; k < 2*NUMWORD-1; k++) {
        for(int i = 0, j = k; i <= k; i++, j--) {
            if (i >= NUMWORD || j >= NUMWORD) continue;

            // uv 계산
            uv = (uint64_t)opa->v[i] * opb->v[j];
            v = (uint32_t)uv;
            u = (uint32_t)(uv >> 32);

            // uv를 r2||r1||r0에 계속 담아주는 작업.
            r0 += v;
            c1 = (r0 < v);
            r1 += u;
            c2 = (r1 < u);
            r1 += c1;
            c2 += (r1 < c1);
            r2 += c2;
        }
        // 한칸 올라감.
        ret->v[k] = r0;
        r0 = r1;
        r1 = r2;
        r2 = 0;
    }
    ret->v[2*NUMWORD-1] = r0;
}

/*  ret = opa mod p , p는 고정이므로 opa의 상위 256비트는 사전계산이 가능하다. 
    사전계산을 잘 정리해서 만든 테이블이 s이다.  
    슈도코드에 맞게 더하거나 뺀 후, 마지막에 캐리 된 만큼 모듈링해준다. 이 때, 반복적으로 p를 빼거나 더하지 말고,
    2p, 3p, 4p, 5p 테이블을 만든 후, 한번에 계산해주면, 연산속도를 높일 수 있다. */
void mod_fast(BN* ret, const BN2* opa)
{
    uint32_t c_add = 0, c_sub = 0;

    // 상위 256비트 사전계산 테이블
    BN s[9] = {
        {opa->v[0] , opa->v[1] , opa->v[2] , opa->v[3] , opa->v[4] , opa->v[5] , opa->v[6] , opa->v[7] },
        {0         , 0         , 0         , opa->v[11], opa->v[12], opa->v[13], opa->v[14], opa->v[15]},
        {0         , 0         , 0         , opa->v[12], opa->v[13], opa->v[14], opa->v[15], 0         },
        {opa->v[8] , opa->v[9] , opa->v[10], 0         , 0         , 0         , opa->v[14], opa->v[15]},
        {opa->v[9] , opa->v[10], opa->v[11], opa->v[13], opa->v[14], opa->v[15], opa->v[13], opa->v[8] },
        {opa->v[11], opa->v[12], opa->v[13], 0         , 0         , 0         , opa->v[8] , opa->v[10]},
        {opa->v[12], opa->v[13], opa->v[14], opa->v[15], 0         , 0         , opa->v[9] , opa->v[11]},
        {opa->v[13], opa->v[14], opa->v[15], opa->v[8] , opa->v[9] , opa->v[10], 0         , opa->v[12]},
        {opa->v[14], opa->v[15], 0         , opa->v[9] , opa->v[10], opa->v[11], 0         , opa->v[13]}};

    // 마지막 계산.
    c_add  = uadd(ret, &s[0], &s[1]);
    c_add += uadd(ret, ret, &s[1]);
    c_add += uadd(ret, ret, &s[2]);
    c_add += uadd(ret, ret, &s[2]);
    c_add += uadd(ret, ret, &s[3]);
    c_add += uadd(ret, ret, &s[4]);
    c_sub =  usub(ret, ret, &s[5]);
    c_sub += usub(ret, ret, &s[6]);
    c_sub += usub(ret, ret, &s[7]);
    c_sub += usub(ret, ret, &s[8]);

    // 모듈링
    if (c_add > c_sub) usub(ret, ret, &nistp256[c_add - c_sub - 1]);
    if (c_add < c_sub) uadd(ret, ret, &nistp256[c_sub - c_add - 1]);
}

/*  EEA의 binary 버전: u와 v를 줄여가며 1 = ax + py 를 만들어가는 것이 목표.
    u = ax1 + py1, v = ax2 + py2    
    u - v = a(x1 - x2) + p(y1 - y2) 따라서 x1 <-- x1 - x2를 한다.  
    u/2 = a(x1/2) + p(y1/2) = a((x1+p)/2) + p((y1+p)/2) since x1 is in Zp   */
void inv(BN* ret, const BN* opa)
{
    BN u = {0, };
    BN v = {0, };
    BN x1 = {0, };
    BN x2 = {0, };
    uint32_t carry = 0;

    set_bn(&u, opa);
    set_bn(&v, &P);
    x1.v[0] = 1;

    while (ucmp(&u, &one) && ucmp(&v, &one)) {
        // u <- u/2 , x1 <- x1/2, 이 때, 홀수이면 LSB가 사라지므로, P를 더해줘서 방지.
        while ((u.v[0] & 1) == 0) {
            rshift1(&u, &u);

            if ((x1.v[0] & 1) == 0) rshift1(&x1, &x1);
            else {                                  
                carry = uadd(&x1, &x1, &P);
                rshift1(&x1, &x1);
                if (carry == 1) x1.v[7] ^= (carry << 31);
            }
        }

        // v <- v/2 , x2 <- x2/2, 이 때, 홀수이면 LSB가 사라지므로, P를 더해줘서 방지.
        while ((v.v[0] & 1) == 0) {
            rshift1(&v, &v);

            if ((x2.v[0] & 1) == 0) rshift1(&x2, &x2);
            else {
                carry = uadd(&x2, &x2, &P);
                rshift1(&x2, &x2);
                if (carry == 1) x2.v[7] ^= (carry << 31);
            }
        }

        // if u > v then u <- u - v and x1 <- x1 - x2 아닌 경우는 반대로
        if (ucmp(&u, &v) >= 0) {
            usub(&u, &u, &v);
            subp(&x1, &x1, &x2);
        } else {
            usub(&v, &v, &u);
            subp(&x2, &x2, &x1);
        }
    }
    
    if (!ucmp(&u, &one)) set_bn(ret, &x1);
    else set_bn(ret, &x2);
}

/*  multiplication in montgomery domain: a x b = (a*b)*R^{-1} mod p
    이 곱셈은 나눗셈이 존재하지 않아 연산속도가 빠르다. 구현방법은 코드의 주석을 참고한다.  */
static void mont(BN* ret, const BN2* opa)
{
    BN2 tmp = {0, };
    BN T_upper = {0, }, T_under = {0, };
    BN U_upper = {0, }, U_under = {0, };

    // T_under || T_upper <-- opa * opb
    memcpy(T_under.v, opa->v, sizeof(BN));
    memcpy(T_upper.v, &opa->v[NUMWORD], sizeof(BN));

    // U <-- T * p' mod R , m' = -p^{-1} --- p는 고정값이므로 사전계산 가능.
    umul_ps(&tmp, &T_under, &P_prime);
    memcpy(&U_under, &tmp, sizeof(BN));
    
    // U_upper || U_under <-- U * p
    umul_ps(&tmp, &U_under, &P);
    memcpy(U_under.v, tmp.v, sizeof(BN));
    memcpy(U_upper.v, &tmp.v[NUMWORD], sizeof(BN));

    // T <-- T + U * p, 이 때 하위 256비트는 무조건 0으로 차고, 버려짐.
    addp(ret, &T_upper, &U_upper);
    if(uadd(&T_under, &T_under, &U_under)) {
        addp(ret, ret, &one);
    }
}

/*  a x b = (ab)R^{-1} mod p
    (a x b) x R^{2} = ((ab)R^{-1}R^2)R^{-1} = ab mod p
    따라서 두번의 몽고메리 곱셈으로 감산이 가능하다.
    이 때, R^2 mod p는 고정값이므로 사전계산을 사용한다. */
void mulp(BN *ret, const BN* opa, const BN* opb)
{
    BN2 T = {0, };

    umul_ps(&T, opa, opb);
    mont(ret, &T);
    umul_ps(&T, ret, &RRmodP);
    mont(ret, &T);
}

/*  montgomery sqr version  */
void sqrp(BN *ret, const BN* opa)
{
    BN2 T = {0, };

    usqr_ps(&T, opa);
    mont(ret, &T);
    umul_ps(&T, ret, &RRmodP);
    mont(ret, &T);
}

    //* https://www.mobilefish.com/services/big_number_equation/big_number_equation.php#equation_output
    //* https://www.boxentriq.com/code-breaking/big-number-calculator



/*

#1. umul_os를 예시로 설명하자면...

opa = a2 || a1 || a0
opb = b2 || b1 || b0

        XX  -- a0 * b0
       XX   -- a1 * b0
   +  XX    -- a2 * b0
   --------
      XXXX  -- r3 || r2 || r1 || r0
   +  
       XX   -- a0 * b1
      XX    -- a1 * b1
   + XX     -- a2 * b1
   --------
     XXXXX  -- r4 || r3 || r2 || r1 || r0 
     
    r2 = 이전 r2 + a0 * b1의 u부분 + a1 * b1의 v부분 -- 이를 알고리즘으로 표현하면,
        uv = ret->v[i+j] + (uint64_t)opa->v[i] * opb->v[j] + u;
        ret->v[i+j] = (uint32_t)uv;
    
    r3에도 이전 u부분을 더해줘야 하므로,
        u = (uint32_t)(uv >> 32);

#2. umul_ps에 대한 설명

위의 umul_os의 그림에서, a1 * b0 와 a0 * b1의 위치가 같다.
이런 식으로 위치가 같은 부분을 전부 더하여 uv에 넣어 계산하는 방식이 ps이다.

ex. r = a0*b0 + 
        a1*b0 + a0*b1 + 
        a2*b0 + a1*b1 + a0*b2
    에서 각 라인을 uv에 담고 전부 더함.
*/