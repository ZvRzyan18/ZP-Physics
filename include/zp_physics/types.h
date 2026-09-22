#ifndef ZP_TYPES_H
#define ZP_TYPES_H

#include <stdint.h>
#include <stddef.h>
#include "zp_physics/hint.h"


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


#ifdef ZP_USE_32BIT_ID
#define ZP_ID_MAX 0xFFFFFFFF
#define ZP_ID_SIZE 4
typedef uint32_t __zp_uid;
#else
#define ZP_ID_MAX 0xFFFF
#define ZP_ID_SIZE 2
typedef uint16_t __zp_uid;
#endif



typedef struct {
 __zp_uid _id;
 uint8_t  _flag;
} zp_handle;

#ifdef zp_has_complex_extension
typedef __complex__ float zp_complex;
#else
typedef zp_vec2 zp_complex;
#endif

typedef zp_vec4 zp_quaternion;

#endif

