#ifndef ZP_VECTOR_H
#define ZP_VECTOR_H

#include "zp_physics/math.h"
#include "zp_physics/types.h"

zp_inline zp_vec2 zp_stv2(const float a) {
 zp_vec2 out;
 out.x = a;
 out.y = a;
 return out;
}

zp_inline zp_vec2 zp_load2(const float *const zp_restrict a) {
 zp_vec2 out;
 out.x = a[0];
 out.y = a[1];
 return out;
}


zp_inline zp_vec2 zp_abs2(const zp_vec2 a) {
 zp_vec2 out;
 out.x = zp_abs(a.x);
 out.y = zp_abs(a.y);
 return out;
}

zp_inline zp_vec2 zp_add2(const zp_vec2 a, const zp_vec2 b) {
 zp_vec2 out;
 out.x = a.x + b.x;
 out.y = a.y + b.y;
 return out;
}

zp_inline zp_vec2 zp_sub2(const zp_vec2 a, const zp_vec2 b) {
 zp_vec2 out;
 out.x = a.x - b.x;
 out.y = a.y - b.y;
 return out;
}


zp_inline zp_vec2 zp_mul2(const zp_vec2 a, const zp_vec2 b) {
 zp_vec2 out;
 out.x = a.x * b.x;
 out.y = a.y * b.y;
 return out;
}


zp_inline zp_vec2 zp_div2(const zp_vec2 a, const zp_vec2 b) {
 zp_vec2 out;
 out.x = a.x / b.x;
 out.y = a.y / b.y;
 return out;
}


zp_inline zp_vec2 zp_fma2(const zp_vec2 a, const zp_vec2 b, const zp_vec2 c) {
 zp_vec2 out;
 out.x = zp_fma(a.x, b.x, c.x);
 out.y = zp_fma(a.y, b.y, c.y);
 return out;
}

zp_pure zp_inline float zp_dot2(const zp_vec2 a, const zp_vec2 b) {
 return zp_fma(a.x, b.x, (a.y * b.y));
}

zp_inline zp_vec2 zp_neg2(const zp_vec2 a) {
 zp_vec2 out;
 out.x = -a.x;
 out.y = -a.y;
 return out;
}

zp_inline zp_vec2 zp_unit2(const zp_vec2 a) {
 return zp_mul2(zp_stv2(zp_rsqrt(zp_dot2(a, a))), a);
}

zp_inline zp_vec2 zp_perp2(const zp_vec2 a) {
 zp_vec2 out;
 out.x = a.y;
 out.y = -a.x;
 return out;
}

zp_inline float zp_cross2(const zp_vec2 a, const zp_vec2 b) {
 return zp_fma(a.x, b.y, -a.y * b.x);
}

zp_inline zp_vec2 zp_cross_sv2(const float a, const zp_vec2 b) {
 zp_vec2 out;
 out.x = -a * b.y;
 out.y = a * b.x;
 return out;
}


zp_inline zp_vec3 zp_stv3(const float a) {
 zp_vec3 out;
 out.x = a;
 out.y = a;
 out.z = a;
 return out;
}

zp_inline zp_vec3 zp_load3(const float *const zp_restrict a) {
 zp_vec3 out;
 out.x = a[0];
 out.y = a[1];
 out.z = a[2];
 return out;
}


zp_inline zp_vec3 zp_abs3(const zp_vec3 a) {
 zp_vec3 out;
 out.x = zp_abs(a.x);
 out.y = zp_abs(a.y);
 out.z = zp_abs(a.z);
 return out;
}


zp_inline zp_vec3 zp_add3(const zp_vec3 a, const zp_vec3 b) {
 zp_vec3 out;
 out.x = a.x + b.x;
 out.y = a.y + b.y;
 out.z = a.z + b.z;
 return out;
}

zp_inline zp_vec3 zp_sub3(const zp_vec3 a, const zp_vec3 b) {
 zp_vec3 out;
 out.x = a.x - b.x;
 out.y = a.y - b.y;
 out.z = a.z - b.z;
 return out;
}


zp_inline zp_vec3 zp_mul3(const zp_vec3 a, const zp_vec3 b) {
 zp_vec3 out;
 out.x = a.x * b.x;
 out.y = a.y * b.y;
 out.z = a.z * b.z;
 return out;
}


zp_inline zp_vec3 zp_div3(const zp_vec3 a, const zp_vec3 b) {
 zp_vec3 out;
 out.x = a.x / b.x;
 out.y = a.y / b.y;
 out.z = a.z / b.z;
 return out;
}


zp_inline zp_vec3 zp_fma3(const zp_vec3 a, const zp_vec3 b, const zp_vec3 c) {
 zp_vec3 out;
 out.x = zp_fma(a.x, b.x, c.x);
 out.y = zp_fma(a.y, b.y, c.y);
 out.z = zp_fma(a.z, b.z, c.z);
 return out;
}

zp_pure zp_inline float zp_dot3(const zp_vec3 a, const zp_vec3 b) {
 return zp_fma(a.x, b.x, zp_fma(a.y, b.y, (a.z * b.z)));
}

zp_inline zp_vec3 zp_neg3(const zp_vec3 a) {
 zp_vec3 out;
 out.x = -a.x;
 out.y = -a.y;
 out.z = -a.z;
 return out;
}

zp_inline zp_vec3 zp_unit3(const zp_vec3 a) {
 return zp_mul3(zp_stv3(zp_rsqrt(zp_dot3(a, a))), a);
}

zp_inline zp_vec3 zp_cross3(const zp_vec3 a, const zp_vec3 b) {
 zp_vec3 out;
 out.x = zp_fma(a.y, b.z, -a.z * b.y);
 out.y = zp_fma(a.z, b.x, -a.x * b.z);
 out.z = zp_fma(a.x, b.y, -a.y * b.x);
 return out;
}


zp_inline zp_vec4 zp_stv4(const float a) {
 zp_vec4 out;
 out.x = a;
 out.y = a;
 out.z = a;
 out.w = a;
 return out;
}

zp_inline zp_vec4 zp_load4(const float *const zp_restrict a) {
 zp_vec4 out;
 out.x = a[0];
 out.y = a[1];
 out.z = a[2];
 out.w = a[3];
 return out;
}

zp_inline zp_vec4 zp_abs4(const zp_vec4 a) {
 zp_vec4 out;
 out.x = zp_abs(a.x);
 out.y = zp_abs(a.y);
 out.z = zp_abs(a.z);
 out.w = zp_abs(a.w);
 return out;
}


zp_inline zp_vec4 zp_add4(const zp_vec4 a, const zp_vec4 b) {
 zp_vec4 out;
 out.x = a.x + b.x;
 out.y = a.y + b.y;
 out.z = a.z + b.z;
 out.w = a.w + b.w;
 return out;
}

zp_inline zp_vec4 zp_sub4(const zp_vec4 a, const zp_vec4 b) {
 zp_vec4 out;
 out.x = a.x - b.x;
 out.y = a.y - b.y;
 out.z = a.z - b.z;
 out.w = a.w - b.w;
 return out;
}


zp_inline zp_vec4 zp_mul4(const zp_vec4 a, const zp_vec4 b) {
 zp_vec4 out;
 out.x = a.x * b.x;
 out.y = a.y * b.y;
 out.z = a.z * b.z;
 out.w = a.w * b.w;
 return out;
}


zp_inline zp_vec4 zp_div4(const zp_vec4 a, const zp_vec4 b) {
 zp_vec4 out;
 out.x = a.x / b.x;
 out.y = a.y / b.y;
 out.z = a.z / b.z;
 out.w = a.w / b.w;
 return out;
}


zp_inline zp_vec4 zp_fma4(const zp_vec4 a, const zp_vec4 b, const zp_vec4 c) {
 zp_vec4 out;
 out.x = zp_fma(a.x, b.x, c.x);
 out.y = zp_fma(a.y, b.y, c.y);
 out.z = zp_fma(a.z, b.z, c.z);
 out.w = zp_fma(a.w, b.w, c.w);
 return out;
}

zp_pure zp_inline float zp_dot4(const zp_vec4 a, const zp_vec4 b) {
 return zp_fma(a.x, b.x, zp_fma(a.y, b.y, zp_fma(a.z, b.z, (a.w * b.w))));
}

zp_inline zp_vec4 zp_neg4(const zp_vec4 a) {
 zp_vec4 out;
 out.x = -a.x;
 out.y = -a.y;
 out.z = -a.z;
 out.w = -a.w;
 return out;
}

zp_inline zp_vec4 zp_unit4(const zp_vec4 a) {
 return zp_mul4(zp_stv4(zp_rsqrt(zp_dot4(a, a))), a);
}



#endif

