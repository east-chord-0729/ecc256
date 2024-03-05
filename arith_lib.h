#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define NUMWORD 8
#define NUMWORD2 16

//* consider :  const BN zero = {0};
//* consider :  constant time implementation
//todo 모든변수는 초기화한다.
//todo ECC 구조체 뭐가 최선인지 생각해보기

/* v: 들어갈 정수: v[n-1] || v[n-2] || ... || v[1] || v[0] 
   c: 캐리
   s: 부호- 양수 = 1, 음수 = -1 */
typedef struct {
    uint32_t v[NUMWORD];
    int32_t s;
} BN;

typedef struct {
    uint32_t v[NUMWORD2];
    int32_t s;
} BN2;

static const BN one  = {{1,0,0,0,0,0,0,0}, 0};
static const BN zero = {{0,0,0,0,0,0,0,0}, 0};

static const BN P = {
    0xffffffff, 0xffffffff, 0xffffffff, 0x00000000, 
    0x00000000, 0x00000000, 0x00000001, 0xffffffff};

static const BN two_inv = {{
    0x00000000, 0x00000000, 0x80000000, 0x00000000, 
    0x00000000, 0x80000000, 0x80000000, 0x7fffffff}, 0};

void set_bn(BN* dest, const BN* src);
int32_t ucmp(const BN* opa, const BN* opb);
uint32_t rshift1(BN* ret, const BN* opa);
uint32_t uadd(BN* ret, const BN* opa, const BN* opb);
uint32_t usub(BN* ret, const BN* opa, const BN* opb);
void umul_os(BN2* ret, const BN* opa, const BN* opb);
void umul_ps(BN2* ret, const BN* opa, const BN* opb);
void usqr_ps(BN2* ret, const BN* opa);

void addp(BN* ret, const BN* opa, const BN* opb);
void subp(BN *ret, const BN* opa, const BN* opb); 
void mod_fast(BN* ret, const BN2* opa);
void mulp(BN *ret, const BN* opa, const BN* opb);
void sqrp(BN *ret, const BN* opa);
void inv(BN* ret, const BN* opa);