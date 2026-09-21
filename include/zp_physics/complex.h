#ifndef ZP_COMPLEX_H
#define ZP_COMPLEX_H

#include "zp_physics/types.h"
#include "zp_physics/math.h"
#include "zp_physics/vector.h"

ZP_CPP_BEGIN

extern zp_complex zp_complex_identity;
extern zp_quaternion zp_quaternion_identity;


/*
 complex for 2d rotation
*/
zp_const zp_inline float zp_creal(const zp_complex c) {
#ifdef zp_has_complex_extension
 return __real__ c;
#else
 return c.x;
#endif
}

zp_const zp_inline float zp_cimag(const zp_complex c) {
#ifdef zp_has_complex_extension
 return __imag__ c;
#else
 return c.y;
#endif
}


zp_inline void zp_set_real(zp_complex *c, float a) {
#ifdef zp_has_complex_extension
 __real__ *c = a;
#else
 c->x = a;
#endif
}

zp_inline void zp_set_imag(zp_complex *c, float a) {
#ifdef zp_has_complex_extension
 __imag__ *c = a;
#else
 c->y = a;
#endif
}



zp_inline zp_complex zp_as_complex2(zp_vec2 a) {
 zp_complex out;
 zp_set_real(&out, a.x);
 zp_set_imag(&out, a.y);
 return out;
}

zp_inline zp_vec2 zp_as_vector2(zp_complex a) {
 zp_vec2 out;
 out.x = zp_creal(a);
 out.y = zp_cimag(a);
 return out;
}


zp_inline zp_complex zp_cabs(zp_complex c) {
#ifdef zp_has_complex_extension
 /*
  return __builtin_cabs(c);
 */
 return zp_as_complex2(zp_abs2(zp_as_vector2(c)));
#else
 return zp_abs2(c);
#endif
}


zp_inline zp_complex zp_cunit(zp_complex c) {
#ifdef zp_has_complex_extension
 return zp_as_complex2(zp_unit2(zp_as_vector2(c)));
#else
 return zp_unit2(c);
#endif
}

zp_const zp_inline zp_complex zp_cmul(const zp_complex a, const zp_complex b) {
#ifdef zp_has_complex_extension
 return a * b;
#else
 zp_complex out;
 out.x = zp_fma(zp_creal(a), zp_creal(b), -(zp_cimag(a) * zp_cimag(b)));
 out.y = zp_fma(zp_creal(a), zp_cimag(b), (zp_cimag(a) * zp_creal(b)));
 return out;
#endif
} 

zp_const zp_inline zp_complex zp_cconj(const zp_complex a) {
#ifdef zp_has_complex_extension
 return __builtin_conjf(a);
#else
 zp_complex out;
 out.x = zp_creal(a);
 out.y = -zp_cimag(a);
 return out;
#endif
} 

zp_complex zp_crotate(const float radians);
float zp_cangle(const zp_complex c);




/*
 quaternion for 3d rotation
*/
zp_const zp_inline zp_quaternion zp_qmul(const zp_quaternion a, const zp_quaternion b) {
 zp_quaternion out;
/*
 out.x = a.x * b.w + a.y * b.z - a.z * b.y + a.w * b.x;
 out.y = -a.x * b.z + a.y * b.w + a.z * b.x + a.w * b.y;
 out.z = a.x * b.y - a.y * b.x + a.z * b.w + a.w * b.z;
 out.w = -a.x * b.x - a.y * b.y - a.z * b.z + a.w * b.w;
*/
 out.x = zp_fma(a.x, b.w, zp_fma(a.y, b.z, -zp_fma(a.z, b.y, a.w * b.x)));
 out.y = zp_fma(-a.x, b.z, zp_fma(a.y, b.w, zp_fma(a.z, b.x, a.w * b.y)));
 out.z = zp_fma(a.x, b.y, -zp_fma(a.y, b.x, zp_fma(a.z, b.w, a.w * b.z)));
 out.w = zp_fma(-a.x, b.x, -zp_fma(a.y, b.y, -zp_fma(a.z, b.z, a.w * b.w)));

 return out;
}


zp_const zp_inline zp_quaternion zp_qconj(const zp_quaternion a) {
 zp_quaternion out;
 out.x = -a.x;
 out.y = -a.y;
 out.z = -a.z;
 out.w = a.w;
 return out;
}

zp_quaternion zp_qfromaxisangle(const float a, const float x, const float y, const float z);
zp_vec4 zp_qtoaxisangle(const zp_quaternion a);


ZP_CPP_END

#endif

