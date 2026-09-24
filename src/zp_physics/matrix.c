#include "zp_physics/matrix.h"
#include "zp_physics/math.h"
#include "zp_physics/complex.h"



zp_mat3x3 zp_mat3x3_identity = {
 .v = {
  {1.0f, 0.0f, 0.0f},
  {0.0f, 1.0f, 0.0f},
  {0.0f, 0.0f, 1.0f}
 }
};

zp_mat4x4 zp_mat4x4_identity = {
 .v = {
  {1.0f, 0.0f, 0.0f, 0.0f},
  {0.0f, 1.0f, 0.0f, 0.0f},
  {0.0f, 0.0f, 1.0f, 0.0f},
  {0.0f, 0.0f, 0.0f, 1.0f}
 }
};


zp_mat3x3 zp_fromquaternion3(const zp_vec4 a) {
 float tmp[9];
 zp_mat3x3 out;

 tmp[0] = a.arr[0] * a.arr[0];
 tmp[1] = a.arr[0] * a.arr[1];
 tmp[2] = a.arr[0] * a.arr[2];
 tmp[3] = a.arr[0] * a.arr[3];

 tmp[4] = a.arr[1] * a.arr[1];
 tmp[5] = a.arr[1] * a.arr[2];
 tmp[6] = a.arr[1] * a.arr[3];

 tmp[7] = a.arr[2] * a.arr[2];
 tmp[8] = a.arr[2] * a.arr[3];

 out.v[0][0] = 1.0f - 2.0f * (tmp[4] + tmp[7]);
 out.v[0][1] =     2.0f * (tmp[1] - tmp[8]);
 out.v[0][2] =     2.0f * (tmp[2] + tmp[6]);

 out.v[1][0] =     2.0f * (tmp[1] + tmp[8]);
 out.v[1][1] = 1.0f - 2.0f * (tmp[0] + tmp[7]);
 out.v[1][2] =     2.0f * (tmp[5] - tmp[3]);

 out.v[2][0] =     2.0f * (tmp[2] - tmp[6]);
 out.v[2][1] =     2.0f * (tmp[5] + tmp[3]);
 out.v[2][2] = 1.0f - 2.0f * (tmp[0] + tmp[4]);

 return out;
}

zp_mat4x4 zp_fromquaternion4(const zp_vec4 a) {
 float tmp[9];
 zp_mat4x4 out;

 tmp[0] = a.arr[0] * a.arr[0];
 tmp[1] = a.arr[0] * a.arr[1];
 tmp[2] = a.arr[0] * a.arr[2];
 tmp[3] = a.arr[0] * a.arr[3];

 tmp[4] = a.arr[1] * a.arr[1];
 tmp[5] = a.arr[1] * a.arr[2];
 tmp[6] = a.arr[1] * a.arr[3];

 tmp[7] = a.arr[2] * a.arr[2];
 tmp[8] = a.arr[2] * a.arr[3];

 out.v[0][0] = 1.0f - 2.0f * (tmp[4] + tmp[7]);
 out.v[0][1] =     2.0f * (tmp[1] - tmp[8]);
 out.v[0][2] =     2.0f * (tmp[2] + tmp[6]);
 out.v[0][3] = 0.0f;
 
 out.v[1][0] =     2.0f * (tmp[1] + tmp[8]);
 out.v[1][1] = 1.0f - 2.0f * (tmp[0] + tmp[7]);
 out.v[1][2] =     2.0f * (tmp[5] - tmp[3]);
 out.v[1][3] = 0.0f;

 out.v[2][0] =     2.0f * (tmp[2] - tmp[6]);
 out.v[2][1] =     2.0f * (tmp[5] + tmp[3]);
 out.v[2][2] = 1.0f - 2.0f * (tmp[0] + tmp[4]);
 out.v[2][3] = 0.0f;

 out.v[3][0] = 0.0f;
 out.v[3][1] = 0.0f;
 out.v[3][2] = 0.0f;
 out.v[3][3] = 1.0f;
 
 return out;
}


/*
 see @complex.c
 
 the idea is based on the Euler's formula 

 e ^ i*θ = cos(θ) + i sin(θ)
 
 e = euler's number (≈2.71828)
 i = imaginary (√-1 or (complex){0, 1})
*/

zp_mat3x3 zp_mat3x3_transform(const zp_complex rotation, const zp_vec2 position) {
 zp_mat3x3 out;
 out.v[0][0] = zp_creal(rotation);
 out.v[0][1] = zp_cimag(rotation);
 out.v[0][2] = 0.0f;
 
 out.v[1][0] = -zp_cimag(rotation);
 out.v[1][1] = zp_creal(rotation);
 out.v[1][2] = 0.0f;
 
 out.v[2][0] = position.x;
 out.v[2][1] = position.y;
 out.v[2][2] = 1.0f;
 return out;
}


zp_mat3x3 zp_transpose3(const zp_mat3x3 a) {
 zp_mat3x3 out;
 out.v[0][0] = a.v[0][0];
 out.v[0][1] = a.v[1][0];
 out.v[0][2] = a.v[2][0];
 
 out.v[1][0] = a.v[0][1];
 out.v[1][1] = a.v[1][1];
 out.v[1][2] = a.v[2][1];

 out.v[2][0] = a.v[0][2];
 out.v[2][1] = a.v[1][2];
 out.v[2][2] = a.v[2][2];
 return out;
}


zp_mat4x4 zp_transpose4(const zp_mat4x4 a) {
 zp_mat4x4 out;
 out.v[0][0] = a.v[0][0];
 out.v[0][1] = a.v[1][0];
 out.v[0][2] = a.v[2][0];
 out.v[0][3] = a.v[3][0];
 
 out.v[1][0] = a.v[0][1];
 out.v[1][1] = a.v[1][1];
 out.v[1][2] = a.v[2][1];
 out.v[1][3] = a.v[3][1];

 out.v[2][0] = a.v[0][2];
 out.v[2][1] = a.v[1][2];
 out.v[2][2] = a.v[2][2];
 out.v[2][3] = a.v[3][2];

 out.v[3][0] = a.v[0][3];
 out.v[3][1] = a.v[1][3];
 out.v[3][2] = a.v[2][3];
 out.v[3][3] = a.v[3][3];
 return out;
}



zp_mat3x3 zp_mulmm3(const zp_mat3x3 a, const zp_mat3x3 b) {
 zp_mat3x3 out;
/*
 out.v[0][0] = a.v[0][0] * b.v[0][0] + a.v[1][0] * b.v[0][1] + a.v[2][0] * b.v[0][2];
 out.v[0][1] = a.v[0][1] * b.v[0][0] + a.v[1][1] * b.v[0][1] + a.v[2][1] * b.v[0][2];
 out.v[0][2] = a.v[0][2] * b.v[0][0] + a.v[1][2] * b.v[0][1] + a.v[2][2] * b.v[0][2];

 out.v[1][0] = a.v[0][0] * b.v[1][0] + a.v[1][0] * b.v[1][1] + a.v[2][0] * b.v[1][2];
 out.v[1][1] = a.v[0][1] * b.v[1][0] + a.v[1][1] * b.v[1][1] + a.v[2][1] * b.v[1][2];
 out.v[1][2] = a.v[0][2] * b.v[1][0] + a.v[1][2] * b.v[1][1] + a.v[2][2] * b.v[1][2];

 out.v[2][0] = a.v[0][0] * b.v[2][0] + a.v[1][0] * b.v[2][1] + a.v[2][0] * b.v[2][2];
 out.v[2][1] = a.v[0][1] * b.v[2][0] + a.v[1][1] * b.v[2][1] + a.v[2][1] * b.v[2][2];
 out.v[2][2] = a.v[0][2] * b.v[2][0] + a.v[1][2] * b.v[2][1] + a.v[2][2] * b.v[2][2];
*/
 out.v[0][0] = zp_fma(a.v[0][0], b.v[0][0], zp_fma(a.v[1][0], b.v[0][1], a.v[2][0] * b.v[0][2]));
 out.v[0][1] = zp_fma(a.v[0][1], b.v[0][0], zp_fma(a.v[1][1], b.v[0][1], a.v[2][1] * b.v[0][2]));
 out.v[0][2] = zp_fma(a.v[0][2], b.v[0][0], zp_fma(a.v[1][2], b.v[0][1], a.v[2][2] * b.v[0][2]));

 out.v[1][0] = zp_fma(a.v[0][0], b.v[1][0], zp_fma(a.v[1][0], b.v[1][1], a.v[2][0] * b.v[1][2]));
 out.v[1][1] = zp_fma(a.v[0][1], b.v[1][0], zp_fma(a.v[1][1], b.v[1][1], a.v[2][1] * b.v[1][2]));
 out.v[1][2] = zp_fma(a.v[0][2], b.v[1][0], zp_fma(a.v[1][2], b.v[1][1], a.v[2][2] * b.v[1][2]));

 out.v[2][0] = zp_fma(a.v[0][0], b.v[2][0], zp_fma(a.v[1][0], b.v[2][1], a.v[2][0] * b.v[2][2]));
 out.v[2][1] = zp_fma(a.v[0][1], b.v[2][0], zp_fma(a.v[1][1], b.v[2][1], a.v[2][1] * b.v[2][2]));
 out.v[2][2] = zp_fma(a.v[0][2], b.v[2][0], zp_fma(a.v[1][2], b.v[2][1], a.v[2][2] * b.v[2][2]));
 return out;
}


zp_mat4x4 zp_mulmm4(const zp_mat4x4 a, const zp_mat4x4 b) {
 zp_mat4x4 out;
/*
 out.v[0][0] = a.v[0][0] * b.v[0][0] + a.v[1][0] * b.v[0][1] + a.v[2][0] * b.v[0][2] + a.v[3][0] * b.v[0][3];
 out.v[0][1] = a.v[0][1] * b.v[0][0] + a.v[1][1] * b.v[0][1] + a.v[2][1] * b.v[0][2] + a.v[3][1] * b.v[0][3];
 out.v[0][2] = a.v[0][2] * b.v[0][0] + a.v[1][2] * b.v[0][1] + a.v[2][2] * b.v[0][2] + a.v[3][2] * b.v[0][3];
 out.v[0][3] = a.v[0][3] * b.v[0][0] + a.v[1][3] * b.v[0][1] + a.v[2][3] * b.v[0][2] + a.v[3][3] * b.v[0][3];

 out.v[1][0] = a.v[0][0] * b.v[1][0] + a.v[1][0] * b.v[1][1] + a.v[2][0] * b.v[1][2] + a.v[3][0] * b.v[1][3];
 out.v[1][1] = a.v[0][1] * b.v[1][0] + a.v[1][1] * b.v[1][1] + a.v[2][1] * b.v[1][2] + a.v[3][1] * b.v[1][3];
 out.v[1][2] = a.v[0][2] * b.v[1][0] + a.v[1][2] * b.v[1][1] + a.v[2][2] * b.v[1][2] + a.v[3][2] * b.v[1][3];
 out.v[1][3] = a.v[0][3] * b.v[1][0] + a.v[1][3] * b.v[1][1] + a.v[2][3] * b.v[1][2] + a.v[3][3] * b.v[1][3];

 out.v[2][0] = a.v[0][0] * b.v[2][0] + a.v[1][0] * b.v[2][1] + a.v[2][0] * b.v[2][2] + a.v[3][0] * b.v[2][3];
 out.v[2][1] = a.v[0][1] * b.v[2][0] + a.v[1][1] * b.v[2][1] + a.v[2][1] * b.v[2][2] + a.v[3][1] * b.v[2][3];
 out.v[2][2] = a.v[0][2] * b.v[2][0] + a.v[1][2] * b.v[2][1] + a.v[2][2] * b.v[2][2] + a.v[3][2] * b.v[2][3];
 out.v[2][3] = a.v[0][3] * b.v[2][0] + a.v[1][3] * b.v[2][1] + a.v[2][3] * b.v[2][2] + a.v[3][3] * b.v[2][3];

 out.v[3][0] = a.v[0][0] * b.v[3][0] + a.v[1][0] * b.v[3][1] + a.v[2][0] * b.v[3][2] + a.v[3][0] * b.v[3][3];
 out.v[3][1] = a.v[0][1] * b.v[3][0] + a.v[1][1] * b.v[3][1] + a.v[2][1] * b.v[3][2] + a.v[3][1] * b.v[3][3];
 out.v[3][2] = a.v[0][2] * b.v[3][0] + a.v[1][2] * b.v[3][1] + a.v[2][2] * b.v[3][2] + a.v[3][2] * b.v[3][3];
 out.v[3][3] = a.v[0][3] * b.v[3][0] + a.v[1][3] * b.v[3][1] + a.v[2][3] * b.v[3][2] + a.v[3][3] * b.v[3][3];
*/
 out.v[0][0] = zp_fma(a.v[0][0], b.v[0][0], zp_fma(a.v[1][0], b.v[0][1], zp_fma(a.v[2][0], b.v[0][2], a.v[3][0] * b.v[0][3])));
 out.v[0][1] = zp_fma(a.v[0][1], b.v[0][0], zp_fma(a.v[1][1], b.v[0][1], zp_fma(a.v[2][1], b.v[0][2], a.v[3][1] * b.v[0][3])));
 out.v[0][2] = zp_fma(a.v[0][2], b.v[0][0], zp_fma(a.v[1][2], b.v[0][1], zp_fma(a.v[2][2], b.v[0][2], a.v[3][2] * b.v[0][3])));
 out.v[0][3] = zp_fma(a.v[0][3], b.v[0][0], zp_fma(a.v[1][3], b.v[0][1], zp_fma(a.v[2][3], b.v[0][2], a.v[3][3] * b.v[0][3])));

 out.v[1][0] = zp_fma(a.v[0][0], b.v[1][0], zp_fma(a.v[1][0], b.v[1][1], zp_fma(a.v[2][0], b.v[1][2], a.v[3][0] * b.v[1][3])));
 out.v[1][1] = zp_fma(a.v[0][1], b.v[1][0], zp_fma(a.v[1][1], b.v[1][1], zp_fma(a.v[2][1], b.v[1][2], a.v[3][1] * b.v[1][3])));
 out.v[1][2] = zp_fma(a.v[0][2], b.v[1][0], zp_fma(a.v[1][2], b.v[1][1], zp_fma(a.v[2][2], b.v[1][2], a.v[3][2] * b.v[1][3])));
 out.v[1][3] = zp_fma(a.v[0][3], b.v[1][0], zp_fma(a.v[1][3], b.v[1][1], zp_fma(a.v[2][3], b.v[1][2], a.v[3][3] * b.v[1][3])));

 out.v[2][0] = zp_fma(a.v[0][0], b.v[2][0], zp_fma(a.v[1][0], b.v[2][1], zp_fma(a.v[2][0], b.v[2][2], a.v[3][0] * b.v[2][3])));
 out.v[2][1] = zp_fma(a.v[0][1], b.v[2][0], zp_fma(a.v[1][1], b.v[2][1], zp_fma(a.v[2][1], b.v[2][2], a.v[3][1] * b.v[2][3])));
 out.v[2][2] = zp_fma(a.v[0][2], b.v[2][0], zp_fma(a.v[1][2], b.v[2][1], zp_fma(a.v[2][2], b.v[2][2], a.v[3][2] * b.v[2][3])));
 out.v[2][3] = zp_fma(a.v[0][3], b.v[2][0], zp_fma(a.v[1][3], b.v[2][1], zp_fma(a.v[2][3], b.v[2][2], a.v[3][3] * b.v[2][3])));

 out.v[3][0] = zp_fma(a.v[0][0], b.v[3][0], zp_fma(a.v[1][0], b.v[3][1], zp_fma(a.v[2][0], b.v[3][2], a.v[3][0] * b.v[3][3])));
 out.v[3][1] = zp_fma(a.v[0][1], b.v[3][0], zp_fma(a.v[1][1], b.v[3][1], zp_fma(a.v[2][1], b.v[3][2], a.v[3][1] * b.v[3][3])));
 out.v[3][2] = zp_fma(a.v[0][2], b.v[3][0], zp_fma(a.v[1][2], b.v[3][1], zp_fma(a.v[2][2], b.v[3][2], a.v[3][2] * b.v[3][3])));
 out.v[3][3] = zp_fma(a.v[0][3], b.v[3][0], zp_fma(a.v[1][3], b.v[3][1], zp_fma(a.v[2][3], b.v[3][2], a.v[3][3] * b.v[3][3])));

 return out;
}



zp_vec3 zp_mulm3v3(const zp_mat3x3 a, const zp_vec3 b) {
 zp_vec3 out;
 out.x = zp_fma(a.v[0][0], b.x, zp_fma(a.v[1][0], b.y, a.v[2][0] * b.z));
 out.y = zp_fma(a.v[0][1], b.x, zp_fma(a.v[1][1], b.y, a.v[2][1] * b.z));
 out.z = zp_fma(a.v[0][2], b.x, zp_fma(a.v[1][2], b.y, a.v[2][2] * b.z));
 return out;
}

zp_vec2 zp_mulm3v2(const zp_mat3x3 a, const zp_vec2 b) {
 zp_vec2 out;
 out.x = zp_fma(a.v[0][0], b.x, zp_fma(a.v[1][0], b.y, a.v[2][0]));
 out.y = zp_fma(a.v[0][1], b.x, zp_fma(a.v[1][1], b.y, a.v[2][1]));
 return out;
}

zp_vec4 zp_mulm4v4(const zp_mat4x4 a, const zp_vec4 b) {
 zp_vec4 out;
 out.x = zp_fma(a.v[0][0], b.x, zp_fma(a.v[1][0], b.y, zp_fma(a.v[2][0], b.z, a.v[3][0] * b.w)));
 out.y = zp_fma(a.v[0][1], b.x, zp_fma(a.v[1][1], b.y, zp_fma(a.v[2][1], b.z, a.v[3][1] * b.w)));
 out.z = zp_fma(a.v[0][2], b.x, zp_fma(a.v[1][2], b.y, zp_fma(a.v[2][2], b.z, a.v[3][2] * b.w)));
 out.w = zp_fma(a.v[0][3], b.x, zp_fma(a.v[1][3], b.y, zp_fma(a.v[2][3], b.z, a.v[3][3] * b.w)));
 return out;
}

zp_vec3 zp_mulm4v3(const zp_mat4x4 a, const zp_vec3 b) {
 zp_vec3 out;
 out.x = zp_fma(a.v[0][0], b.x, zp_fma(a.v[1][0], b.y, zp_fma(a.v[2][0], b.z, a.v[3][0])));
 out.y = zp_fma(a.v[0][1], b.x, zp_fma(a.v[1][1], b.y, zp_fma(a.v[2][1], b.z, a.v[3][1])));
 out.z = zp_fma(a.v[0][2], b.x, zp_fma(a.v[1][2], b.y, zp_fma(a.v[2][2], b.z, a.v[3][2])));
 return out;
}


