#ifndef ZP_MATH_H
#define ZP_MATH_H

#include "zp_physics/hint.h"
#include <math.h>

ZP_CPP_BEGIN

typedef union {
	float f;
 uint32_t i;
 int32_t s;
} __f_floatbits;


zp_const zp_inline float zp_nan(void) {
 __f_floatbits bits;
 bits.i = 0x7F800001;
 return bits.f;
}



zp_const zp_inline float zp_inf(void) {
 __f_floatbits bits;
 bits.i = 0x7F800000;
 return bits.f;
}

/*
 not usually a fastest approach, but still usefull
*/
zp_const zp_inline float zp_copysign(float a, float b) {
 __f_floatbits a_bits, b_bits;
 a_bits.f = a;
 b_bits.f = b;
 a_bits.i = (a_bits.i & 0x7FFFFFFF) | (b_bits.i & 0x80000000); 
 return a_bits.f;
}

/*
 use std fmaf instead?
*/
zp_const zp_inline float zp_fma(const float a, const float b, const float c) {
/*
#if defined(__aarch64__)
 float result;
 __asm__ volatile(
 "fmadd %s0, %s1, %s2, %s3"
 : "=w"(result)
 : "w"(a), "w"(b), "w"(c)
 );
 return result;
#else
*/
 return fmaf(a, b, c);
}

/*
 use std fabsf instead?
*/
zp_const zp_inline float zp_abs(const float x) {
#if defined(__aarch64__)
 float result;
 __asm__ volatile(
 "fabs %s0, %s1"
 : "=w"(result)
 : "w"(x)
 );
 return result;
#else
 __f_floatbits f;
 f.f = x;
 f.i &= 0x7FFFFFFF;
 return f.f;
#endif
}


zp_const zp_inline float zp_min(const float a, const float b) {
#if defined(__aarch64__)
 float result;
 __asm__ volatile(
 "fmin %s0, %s1, %s2"
 : "=w"(result)
 : "w"(a), "w"(b)
 );
 return result;
#else
 return (a < b) ? a : b;
#endif
}

zp_const zp_inline float zp_max(const float a, const float b) {
#if defined(__aarch64__)
 float result;
 __asm__ volatile(
 "fmax %s0, %s1, %s2"
 : "=w"(result)
 : "w"(a), "w"(b)
 );
 return result;
#else
 return (a > b) ? a : b;
#endif
}

zp_const uint64_t zp_alignto(const uint64_t x, const uint8_t alignment);
zp_const uint32_t zp_nextp2(const uint32_t x);

/*
 use std math instead? 
 for round, ceil, floor
*/
zp_pure float zp_round(const float x);
zp_pure float zp_ceil(const float x);
zp_pure float zp_floor(const float x);

/*
 trigs
*/
void zp_sincos(const float theta, float *const zp_restrict s, float *const zp_restrict c);
zp_pure float zp_sin(const float theta);
zp_pure float zp_cos(const float theta);
float zp_tan(const float theta);

/*
 inv trig
*/
zp_pure float zp_asin(const float x);
zp_pure float zp_acos(const float x);
zp_pure float zp_atan(const float x);
zp_pure float zp_atan2(const float y, const float x);

/*
 exponential and logarithm
*/
zp_pure float zp_log2(const float x);
zp_pure float zp_exp2(const float x);
zp_pure float zp_tanh(const float x);
zp_pure float zp_rsqrt(const float x);
zp_pure float zp_sqrt(const float x);


ZP_CPP_END

#endif

