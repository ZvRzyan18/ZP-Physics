#ifndef ZP_TYPES_H
#define ZP_TYPES_H

#include <stdint.h>
#include <stddef.h>
#include "zp_physics/hint.h"

/*
 NOTE : DONT OVERLAP THE VALUES POINTER, DOING SO MIGHT RESULT IN UNDEFINED BEHAVIOUR
*/
typedef union {
 struct {	float x, y; };
 float arr[2];
} zp_vec2;

typedef union {
 struct {	float x, y, z; };
 float arr[3];
} zp_vec3;

typedef union {
 struct {	float x, y, z, w; };
 float arr[4];
} zp_vec4;

typedef union {
 zp_vec3 basis[3];
 float v[3][3];
 float arr[9];
} zp_mat3x3;

typedef union {
 zp_vec4 basis[4];
 float v[4][4];
 float arr[16];
} zp_mat4x4;



typedef uint32_t zp_handle;

#ifdef zp_has_complex_extension
typedef __complex__ float zp_complex;
#else
typedef zp_vec2 zp_complex;
#endif

typedef zp_vec4 zp_quaternion;

#endif

