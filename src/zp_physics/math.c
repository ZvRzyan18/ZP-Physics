#include "zp_physics/math.h"

/*

 Deterministic math. i guess?

*/


/*
 next multiple of
*/
zp_const uint64_t zp_alignto(const uint64_t x, const uint8_t alignment) {
 return (x + (uint64_t)alignment - 1) & ~((uint64_t)alignment - 1);
}

/*
 next power of 2
*/
zp_const uint32_t zp_nextp2(const uint32_t x) {
 uint32_t mx = x;
 mx--;
 mx |= mx >> 1;
 mx |= mx >> 2;
 mx |= mx >> 4;
 mx |= mx >> 8;
 mx |= mx >> 16;
 return ++mx;
}

/*
 2^23, large enough to completly round off fractionals
*/
static const float HUGE_NUM[4] = {
 8388608.0f, -8388608.0f,
 -8388608.0f, 8388608.0f,
};
/*
 round to nearest, every .5 it increases.
*/
zp_pure float zp_round(const float x) {
	__f_floatbits bits;
	bits.f = x;
	const float*const huge = HUGE_NUM + ((bits.i >> 31) << 1);
	return (x + huge[0]) + huge[1];
}

/*
 every fractional digits, automatically rounds up towards positive infinity
*/
zp_pure float zp_ceil(const float x) {
 int i = (int)x;
 return (float)(x > 0 && x != (float)i) ? i + 1 : i;
}

/*
 every fractional digits, automatically rounds up towards negative infinity
*/
zp_pure float zp_floor(const float x) {
 int i = (int)x;
 return (float)(x < 0 && x != (float)i) ? i - 1 : i;
}


/*
 coeffs polynomial with lowest degree
 just acceptable enough for rough estimation
*/
/* original coeffs */
/*
  3.515696e-02f,
 -4.909662e-01f,
  1.00000000e-00f,

static const float COSINE[3] = {

//  scaled coeffs
 2.97037352e-18f,
-4.51285593e-09f,
 1.00000000e-00f,
};
*/

/*
 scaled by the 
 COSINE[0] *= powf(M_PI*2 / 65536.0, 4)
 COSINE[1] *= powf(M_PI*2 / 65536.0, 2)

*/

/*-------------------------------------------------------------------*/

/*

               point in unit circle
               *
              / |
             /  |
  radius -> /   |
           /    | -> sin
          /     |
         / |    |
 origin *------ *
       angle  | -> cos

  cos(theta) = x-axis of theta angle
  sin(theta) = y-axis of theta angle
*/ 
 /*
  convert to fixed point of (65536 / TAU)
  and use intentional overflow so we dont have to
  use remainder operation.
 */
 /*
void zp_sincos(const float theta, float *const zp_restrict s, float *const zp_restrict c) {

 uint16_t tx = (uint16_t)(uint32_t)(zp_abs(theta) * 10430.378350f);
 int16_t ix = (int16_t)tx;
 uint8_t fp = (tx >= 0x3FFF) && (tx <= 0xBFFF);
 ix = fp ? 0x7FFF - ix : ix;
 // convert it back
 // np_real t = ((np_real)ix) * 0.0000958738f;
 
 float t = ((float)ix);

 float x2 = t * t;
 float mx = zp_fma(zp_fma(COSINE[0], x2, COSINE[1]), x2, COSINE[2]);
 float mx1 = zp_sqrt(1.0f - mx * mx);
 *c = fp ? -mx : mx;
 *s = ((theta < 0.0f) ^ (tx > 0x7FFF)) ? -mx1 : mx1;
}

zp_pure float zp_sin(const float theta) { 
 uint16_t tx = 0x3FFF - (uint16_t)(uint32_t)(theta * 10430.378350f);
 int16_t ix = tx;
 uint8_t fp = (tx >= 0x3FFF) && (tx <= 0xBFFF);
 ix = fp ? 0x7FFF - ix : ix;
 float t = ((float)ix);
 float x2 = t * t;
 float mx = zp_fma(zp_fma(COSINE[0], x2, COSINE[1]), x2, COSINE[2]);
 return fp ? -mx : mx;
}



zp_pure float zp_cos(const float theta) { 
 uint16_t tx = (uint16_t)(uint32_t)(theta * 10430.378350f);
 int16_t ix = tx;
 uint8_t fp = (tx >= 0x3FFF) && (tx <= 0xBFFF);
 ix = fp ? 0x7FFF - ix : ix;
 float t = ((float)ix);
 float x2 = t * t;
 float mx = zp_fma(zp_fma(COSINE[0], x2, COSINE[1]), x2, COSINE[2]);
 return fp ? -mx : mx;
}
*/


/*
 approximation 1 degree
 interval : [0.0:0.785398163397448309615660845820]
 p(x) = (cos(sqrt(x))-1) / x
 cos(x) = 1 + x2 * p(x2)
 
 coeffs
  0.040591002033280692f,
 -0.49989513679207137f,
 
 scaled 
 0.040591002033280692f * powf(M_PI*2 / 65536.0, 4)
 -0.49989513679207137f * powf(M_PI*2 / 65536.0, 2)
 
*/


/*
 new version of trig function.
 slightly more expensive but it encreases the accuracy
*/
static float COSINE[2] = {
  3.429489862425e-18f,
 -4.594928810527e-09f,
};

void zp_sincos(const float theta, float *const zp_restrict s, float *const zp_restrict c) {
 uint16_t tx = (uint16_t)(uint32_t)(zp_abs(theta) * 10430.378350f);
 
 int16_t r = (int16_t)tx;
 uint8_t fp = (tx >= 0x3FFF) && (tx <= 0xBFFF);
 r = fp ? 0x7FFF - r : r;
 uint8_t flip = r > 8192;
 r = flip ? (16384 - r) : r;
 
 float mx = ((float)r);
 float x2 = mx * mx;
 float cosine = zp_fma(x2, zp_fma(COSINE[0], x2, COSINE[1]), 1.0f);
 float sine = zp_sqrt(1.0f - cosine * cosine);
 
 float os, oc;
 os = flip ? cosine : sine;
 oc = flip ? sine : cosine;
 os = ((theta < 0.0f) ^ (tx > 0x7FFF)) ? -os : os;
 oc = fp ? -oc : oc;
 *s = os;
 *c = oc;
}



zp_pure float zp_sin(const float theta) { 
 uint16_t tx = (uint16_t)(uint32_t)(zp_abs(theta) * 10430.378350f);
 
 int16_t r = (int16_t)tx;
 uint8_t fp = (tx >= 0x3FFF) && (tx <= 0xBFFF);
 r = fp ? 0x7FFF - r : r;
 uint8_t flip = r > 8192;
 r = flip ? (16384 - r) : r;
 
 float mx = ((float)r);
 float x2 = mx * mx;
 float cosine = zp_fma(x2, zp_fma(COSINE[0], x2, COSINE[1]), 1.0f);
 float sine = zp_sqrt(1.0f - cosine * cosine);
 
 float os;
 os = flip ? cosine : sine;
 return ((theta < 0.0f) ^ (tx > 0x7FFF)) ? -os : os;
}



zp_pure float zp_cos(const float theta) { 
 uint16_t tx = (uint16_t)(uint32_t)(zp_abs(theta) * 10430.378350f);
 
 int16_t r = (int16_t)tx;
 uint8_t fp = (tx >= 0x3FFF) && (tx <= 0xBFFF);
 r = fp ? 0x7FFF - r : r;
 uint8_t flip = r > 8192;
 r = flip ? (16384 - r) : r;
 
 float mx = ((float)r);
 float x2 = mx * mx;
 float cosine = zp_fma(x2, zp_fma(COSINE[0], x2, COSINE[1]), 1.0f);
 float sine = zp_sqrt(1.0f - cosine * cosine);
 
 float oc;
 oc = flip ? sine : cosine;
 return fp ? -oc : oc;
}



float zp_tan(const float theta) {
	float s, c;
 zp_sincos(theta, &s, &c);
 return s / c;
}


/*
 arc sine
          • (x, y) = normalized value
         /|
        / |
       /  |
      /   |
     /__  |
    /   | |
   •-------
   ^
  angle
 
  angle = asin(y)
*/

/*
 approximation

 acos(x) = acos(x)/sqrt(1.0 - x), interval : [0, 1]
 asin(x) = pi/2 - acos(x)
*/
static const float ASIN[3] = {
  0.46707077e-01f,
 -0.20212058e-00f,
  1.57021169e-00f,
};

zp_pure float zp_asin(const float x) {
 float mx, mx1;
 mx = zp_min(zp_abs(x), 1.0f);
 mx1 = zp_fma(ASIN[0], mx, ASIN[1]);
 mx1 = zp_fma(mx1, mx, ASIN[2]);
 mx1 = mx1 * zp_sqrt(1.0f - mx);
 return (1.570796f - ((x < 0.0f) ? 3.14159f - mx1 : mx1));
}

/*
 arc cosine
          • (x, y) = normalized value
         /|
        / |
       /  |
      /   |
     /__  |
    /   | |
   •-------
   ^
  angle
 
  angle = acos(x)
*/
/*
 approxiamtion

 acos(x) = acos(x)/sqrt(1.0 - x), interval : [0, 1]
*/

static const float ACOS[2] = {
 -0.15658276e-00f,
  1.56496682e-00f,
};

zp_pure float zp_acos(const float x) {
 float mx, mx1;
 mx = zp_min(zp_abs(x), 1.0f);
 mx1 = zp_fma(ACOS[0], mx, ACOS[1]);
 mx1 = mx1 * zp_sqrt(1.0f - mx);
 return ((x < 0.0f) ? 3.14159f - mx1 : mx1);
}

/*
 arc tangent
          • (x, y) normalized value
         /|
        / |
       /  |
      /   |
     /__  |
    /   | |
   •-------
     ^
    angle
 
  angle = atan(y/x)
*/

static const float ATAN[2] = {
/*
 -6.038020e-02f,
  1.770501e-01f,
 -3.323025e-01f,
*/
  1.187314e-01f,
 -3.257921e-01f
};

/*
 # arc tangent
                     ________
 atan(x) = asin(x / √1 + x*x ) 

 • alternative (only works on positive x):
                     ________
 atan(x) = acos(1 / √1 + x*x )
*/

/*
 implementation
 
 for |x| <= 1.0
 
 atanPoly(x) = c0 * x + c1 * x² + c2 * x³....
 
 for |x| >= 1.0
 
 atan(x) pi_half - atanPoly(1.0 / x)
*/

zp_pure float zp_atan(const float x) {
 float x2, mx, out;
 uint8_t hi;
 
 mx = zp_abs(x);
 hi = mx > 1.0f;
 if(!hi) {
  x2 = mx * mx;
  out = zp_fma((mx * x2), zp_fma(ATAN[0], x2, ATAN[1]), mx);
  return zp_copysign(mx, x);
 }
 mx = 1.0f / mx;
 x2 = mx * mx;
 out = zp_fma((mx * x2), zp_fma(ATAN[0], x2, ATAN[1]), mx);
 out = 1.570796f - out;
 return zp_copysign(out, x);
}

/*
 arc tangent but with separate sine and cosine input and quadrant correction
 because sin/cos = tan
*/
zp_pure float zp_atan2(const float y, const float x) {
 float ratio, mx, x2;
 float ax = zp_abs(x);
 float ay = zp_abs(y);

 if(ax >= ay) {
  ratio = ay / ax;
  mx = ratio;
  x2 = mx * mx;
  mx = zp_fma((mx * x2), zp_fma(ATAN[0], x2, ATAN[1]), mx);
 } else {
  ratio = ax / ay;
  mx = ratio;
  x2 = mx * mx;
  mx = 1.570796f - zp_fma((mx * x2), zp_fma(ATAN[0], x2, ATAN[1]), mx);
 }
    
 mx = (x < 0.0f) ? (3.14159f - mx) : mx;
 mx = zp_copysign(mx, y);
 return mx;
}


/*-------------------------------------------------------------------*/
/*
 log base 2 coeffs
 accurate enough for at least 2-3 decimal digits.
*/
static const float LOG2[5] = {
 -8.161580e-02f,
  6.451423e-01f,
 -2.120675e-00f,
  4.070090e-00f,
 -2.512854e-00f,
};
/*

 
                                        ______
                                    __ -
  logarithmic growth ->         __ -
                            _ -
                           /
                          -
                         /
    udefined negatives  |
      -------|-------|--}* -> point exactly to one
                       /
                       | -> below 1.0 are negative values
                       /
   
-> gradually goes to zero with negative sign


  log2(x) = ln(x) / ln(2)
 
*/
zp_pure float zp_log2(const float x) {
 __f_floatbits mantissa, bits;
 float ma, mx;
 
 if(x >= 1.0f && x < 2.0f) {
  ma = x;
  return ((((LOG2[0] * ma + LOG2[1]) * ma + LOG2[2]) * ma + LOG2[3]) * ma + LOG2[4]);
 } else if(x < 1.0f) {
 	bits.f = 1.0f / x;
  mantissa.i = 1065353216U | (bits.i & 0x007FFFFF);
  ma = mantissa.f;
  mx = ((((LOG2[0] * ma + LOG2[1]) * ma + LOG2[2]) * ma + LOG2[3]) * ma + LOG2[4]);
  mx = (float)((int32_t)(bits.i >> 23) - 127) + mx;
  return -mx;
 }
 
 bits.f = x;
 mantissa.i = 1065353216U | (bits.i & 0x007FFFFF);
 ma = mantissa.f;
 mx = ((((LOG2[0] * ma + LOG2[1]) * ma + LOG2[2]) * ma + LOG2[3]) * ma + LOG2[4]);
 mx = (float)((int32_t)(bits.i >> 23) - 127) + mx;
 return mx;
}

/*-------------------------------------------------------------------*/
/*
 exponential base 2 coeffs
 accurate enough for at least 2-3 decimal digits.
*/
static const float EXP2[5] = {
 0.96857114e-02f,
 0.55922025e-01f,
 0.24021865e-00f,
 0.69312103e-00f,
 1.00000015e-00f,
};

/*

 
                exponential growth

                        /
                       |
                      /   
                     |
                    /
                ---*  <- point exactly to one
              __
    _______---      <- negative values are always below one

    exp2(x) = exp(x * ln(2))
*/

zp_pure float zp_exp2(const float x) { 
 float mx, a;
	__f_floatbits bits;
 int32_t whole;

 if(zp_abs(x) < 0.5f) {
 	a = x;
  return ((((EXP2[0] * a + EXP2[1]) * a + EXP2[2]) * a + EXP2[3]) * a + EXP2[4]);
 }
	a = x;
 whole = (float)(a - 0.5f);
 a = a - ((float)whole);
 bits.i = (uint32_t)((int32_t)(127 + whole) << 23);
 mx = ((((EXP2[0] * a + EXP2[1]) * a + EXP2[2]) * a + EXP2[3]) * a + EXP2[4]);
 return bits.f * mx;
}

/*-------------------------------------------------------------------*/

/*
 hyperbolic tangent 

 1.0 ---------------------------------
                    ____
                ___
            ___
         __
       __
     _
   -
  -
-------------------------------------
 for negative values it is just '-tanh(abs(x))'
 
 as the value increases it slowly approaches 1.0
 useful for smooth transition with clamped limits.

 tanh(x) = (e ^ x - 1.0) / (e ^ x + 1.0);

*/
/*
zp_pure float zp_tanh(const float x) {
 if(zp_abs(x) >= 45.0f)
  return zp_copysign(1.0f, x);
 const float log2_e = 1.44269504f;
 float exponential = zp_exp2(2.0 * x * log2_e);
 return (exponential - 1.0) / (exponential + 1.0);
}
*/
static const float TANH[6] = {
 9.49274892e-02f,
-3.29312635e-01f,

 4.78187872e-02f,
-3.72573327e-01f,
 9.85379375e-01f,
 1.03154878e-01f,
};

zp_pure float zp_tanh(const float x) {
 float mx = zp_abs(x);
 if(mx < 1.0) {
  float x2 = mx * mx;
  float mm = zp_fma((x2 * mx), zp_fma(TANH[0], x2, TANH[1]), mx);
  return zp_copysign(mm, x);
 } else if(mx < 3.0) {
  float mm = zp_fma(zp_fma(zp_fma(TANH[2], mx, TANH[3]), mx, TANH[4]), mx, TANH[5]);
  return zp_copysign(mm, x);
 } 
 return zp_copysign(1.0f, x);
}



#if defined(__aarch64__)

zp_pure float zp_rsqrt(const float x) {
 float out;
 __asm__ volatile (
  "fsqrt %s0, %s1\n"
  : "=w" (out)
  : "w" (x)
 );
 return 1.0f / out;
}

zp_pure float zp_sqrt(const float x) {
 float out;
 __asm__ volatile (
  "fsqrt %s0, %s1\n"
  : "=w" (out)
  : "w" (x)
 );
 return out;
}

#elif defined(__x86_64__)

zp_pure float zp_rsqrt(const float x) {
 float o;
 __asm__ volatile (
  "movss %[in], %%xmm0 \n"
  "sqrtss %%xmm0, %%xmm0 \n"
  "movss %%xmm0, %[out]"
  : [out] "=m" (o)
  : [in] "m" (x)
  : "xmm0"
 );
 return 1.0f / o;
}

zp_pure float zp_sqrt(const float x) {
 float o;
 __asm__ volatile (
  "movss %[in], %%xmm0 \n"
  "sqrtss %%xmm0, %%xmm0 \n"
  "movss %%xmm0, %[out]"
  : [out] "=m" (o)
  : [in] "m" (x)
  : "xmm0"
 );
 return o;
}

#else

/*
 Quake III fast inv sqrt
*/
zp_pure float zp_rsqrt(const float x) {
 float f, mx, x_half;
 __f_floatbits bits, x_bits;
 x_bits.f = x;
 bits.i =   0x5F378171 - (x_bits.i >> 1);
 mx =       bits.f;
 x_half =   x * 0.5f;
 f =        mx;
 f =        (f * (1.5f - (x_half * f * f)));
 f =        (f * (1.5f - (x_half * f * f)));
 return f;
}


zp_pure float zp_sqrt(const float x) {
 float f, mx, x_half;
 __f_floatbits bits, x_bits;
 x_bits.f = x;
 bits.i =   0x5F378171 - (x_bits.i >> 1);
 mx =       bits.f;
 x_half =   x * 0.5f;
 f =        mx;
 f =        (f * (1.5f - (x_half * f * f)));
 f =        (f * (1.5f - (x_half * f * f)));
 return x * f;
}

#endif


