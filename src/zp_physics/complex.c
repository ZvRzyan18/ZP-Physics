#include "zp_physics/complex.h"

zp_complex zp_complex_identity = {.x = 1.0f, .y = 0.0f};
zp_quaternion zp_quaternion_identity = {.x = 0.0f, .y = 0.0f, .z = 0.0f, .w = 1.0f};


zp_complex zp_crotate(const float radians) {
 zp_complex out;
 zp_sincos(radians, &out.y, &out.x);
 return out;
}


float zp_cangle(const zp_complex c) {
 return zp_atan2(c.y, c.x);
}



zp_quaternion zp_qfromaxisangle(const float a, const float x, const float y, const float z) {
 zp_quaternion out;
 float s;
 zp_sincos(a * 0.5f, &s, &out.w);
	out.x = x * s;
 out.y = y * s;
	out.z = z * s;
	return out;
}

zp_vec4 zp_qtoaxisangle(const zp_quaternion a) {
 zp_vec4 out;
 out.w = zp_asin(a.w);
 out.w += out.w;
 const float b = zp_rsqrt(1.0f - a.w * a.w);
 out.x = a.x * b;
 out.y = a.y * b;
 out.z = a.z * b;
 return out;
}




