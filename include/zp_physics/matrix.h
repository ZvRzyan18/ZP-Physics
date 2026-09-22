#ifndef ZP_MATRIX_H
#define ZP_MATRIX_H

#include "zp_physics/types.h"
#include "zp_physics/hint.h"


#ifdef __cplusplus
extern "C" {
#endif


extern zp_mat3x3 zp_mat3x3_identity;
extern zp_mat4x4 zp_mat4x4_identity;

zp_mat3x3 zp_fromquaternion3(const zp_vec4 a);
zp_mat4x4 zp_fromquaternion4(const zp_vec4 a);

zp_mat3x3 zp_mat3x3_transform(const zp_complex rotation, const zp_vec2 position);

zp_mat3x3 zp_transpose3(const zp_mat3x3 a);
zp_mat4x4 zp_transpose4(const zp_mat4x4 a);

zp_mat3x3 zp_mulmm3(const zp_mat3x3 a, const zp_mat3x3 b);
zp_mat4x4 zp_mulmm4(const zp_mat4x4 a, const zp_mat4x4 b);

zp_vec3 zp_mulm3v3(const zp_mat3x3 a, const zp_vec3 b);
zp_vec2 zp_mulm3v2(const zp_mat3x3 a, const zp_vec2 b);
zp_vec4 zp_mulm4v4(const zp_mat4x4 a, const zp_vec4 b);
zp_vec3 zp_mulm4v3(const zp_mat4x4 a, const zp_vec3 b);



#ifdef __cplusplus
}
#endif



#endif



