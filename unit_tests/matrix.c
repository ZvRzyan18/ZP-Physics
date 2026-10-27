#include "matrix.h"
#include "oxb1tz/matrix.h"
#include "oxb1tz/vector.h"


#ifdef OXB_TEST

oxb_inline float MATRIX_TEST_RANDOM(float min, float max) {
 return min + ((float)rand() / (float)RAND_MAX) * (max - min);
}

static void TEST_MAT3X3_INV(float m[9], float inv[9]) {
 float a00 = m[0], a01 = m[3], a02 = m[6];
 float a10 = m[1], a11 = m[4], a12 = m[7];
 float a20 = m[2], a21 = m[5], a22 = m[8];

 float c00 = a11 * a22 - a12 * a21;
 float c01 = a02 * a21 - a01 * a22;
 float c02 = a01 * a12 - a02 * a11;
 float det = a00 * c00 + a10 * c01 + a20 * c02;
 float invDet = 1.0f / det;
 inv[0] = c00 * invDet;
 inv[1] = c01 * invDet;
 inv[2] = c02 * invDet;

 inv[3] = (a12 * a20 - a10 * a22) * invDet;
 inv[4] = (a00 * a22 - a02 * a20) * invDet;
 inv[5] = (a02 * a10 - a00 * a12) * invDet;

 inv[6] = (a10 * a21 - a11 * a20) * invDet;
 inv[7] = (a01 * a20 - a00 * a21) * invDet;
 inv[8] = (a00 * a11 - a01 * a10) * invDet;
}

static void MAT3_MULV(float out[3], float a[3], float b[9]) {
 for(int row = 0; row < 3; row++) {
  out[row] = 0.0f;
  for(int col = 0; col < 3; col++) {
   out[row] += a[col] * b[col * 3 + row];
  }
 }
}

static void mat3x3_test() {
 const float epsilon = 1e-2f;
 const float min_range = -100.0f;
 const float max_range =  100.0f;
 const int loop_repeat = 100;
 
 printf("--transpose test\n");
 /* transpose test */
 for(int kk = 0; kk < loop_repeat; kk++) {
  float tmp[9];
  
  for(int i = 0; i < 9; i++) {
  	tmp[i] = MATRIX_TEST_RANDOM(min_range, max_range);
  }
  
  oxb_mat3x3 aa;
  aa._basis[0] = oxb_vec3_load(tmp);
  aa._basis[1] = oxb_vec3_load(tmp+3);
  aa._basis[2] = oxb_vec3_load(tmp+6);
  aa = oxb_mat3x3_transpose(aa);
  
  float transposed[9];
  float transposed1[9];
  oxb_mat3x3_tofloat(transposed1, aa);
  for(int i = 0; i < 3; i++) {
   for(int j = 0; j < 3; j++) {
    transposed[j * 3 + i] = tmp[i * 3 + j];
   }
  }
  
  for(int i = 0; i < 9; i++) {
  /*
   if(fabsf(transposed1[i] - transposed[i]) < epsilon)
    printf("%f, %f\n", transposed1[i], transposed[i]);
   */
   assert(fabsf(transposed1[i] - transposed[i]) < epsilon);
  }
 }
 
 printf("--mul test\n");
 /* mul test */
 for(int kk = 0; kk < loop_repeat; kk++) {
  float tmp[9], tmp1[9];
  
  for(int i = 0; i < 9; i++) {
  	tmp[i] = MATRIX_TEST_RANDOM(min_range, max_range);
  	tmp1[i] = MATRIX_TEST_RANDOM(min_range, max_range);
  }

  oxb_mat3x3 aa;
  aa._basis[0] = oxb_vec3_load(tmp);
  aa._basis[1] = oxb_vec3_load(tmp+3);
  aa._basis[2] = oxb_vec3_load(tmp+6);

  oxb_mat3x3 bb;
  bb._basis[0] = oxb_vec3_load(tmp1);
  bb._basis[1] = oxb_vec3_load(tmp1+3);
  bb._basis[2] = oxb_vec3_load(tmp1+6);
  
  oxb_mat3x3 rr = oxb_mat3x3_mul(aa, bb);
  float product1[9], product2[9];
  oxb_mat3x3_tofloat(product1, rr);

  for(int col = 0; col < 3; col++) {
   for(int row = 0; row < 3; row++) {
    product2[row + col * 3] = 0.0f;
    for(int k = 0; k < 3; k++) {
     product2[row + col * 3] = fmaf(tmp[row + k * 3], tmp1[k + col * 3], product2[row + col * 3]);
    }
   }
  }

  for(int i = 0; i < 9; i++) {
   assert(fabsf(product1[i] - product2[i]) < epsilon);
  }
 }
 
 printf("--inv test\n");
 /* inv test */
 for(int kk = 0; kk < loop_repeat; kk++) {
  float tmp[9];
  
  for(int i = 0; i < 9; i++) {
  	tmp[i] = MATRIX_TEST_RANDOM(min_range, max_range);
  }
  
  oxb_mat3x3 aa;
  aa._basis[0] = oxb_vec3_load(tmp);
  aa._basis[1] = oxb_vec3_load(tmp+3);
  aa._basis[2] = oxb_vec3_load(tmp+6);
  aa = oxb_mat3x3_inv(aa);
  
  float r1[9];
  float r2[9];
  oxb_mat3x3_tofloat(r1, aa);
  TEST_MAT3X3_INV(tmp, r2);
  
  for(int i = 0; i < 9; i++) {
   if((isinf(r1[i]) || isnan(r1[i])) && (isinf(r2[i]) || isnan(r2[i])))
    continue;
   assert(fabsf(r1[i] - r2[i]) < epsilon);
  }
 }
 
 printf("--vector matrix mul test\n");
 /* inv test */
 for(int kk = 0; kk < loop_repeat; kk++) {
  float tmp[9];
  float tmp1[3];
  for(int i = 0; i < 3; i++) {
  	tmp1[i] = MATRIX_TEST_RANDOM(min_range, max_range);
  }
  for(int i = 0; i < 9; i++) {
  	tmp[i] = MATRIX_TEST_RANDOM(min_range, max_range);
  }
  
  oxb_mat3x3 aa;
  aa._basis[0] = oxb_vec3_load(tmp);
  aa._basis[1] = oxb_vec3_load(tmp+3);
  aa._basis[2] = oxb_vec3_load(tmp+6);
  oxb_vec3 bb = oxb_vec3_load(tmp1);
  bb = oxb_mat3x3_mulvec(bb, aa);
  
  float cc[3], dd[3];
  MAT3_MULV(cc, tmp1, tmp);
  oxb_vec3_store(dd, bb);
  for(int i = 0; i < 3; i++) {
   if((isinf(dd[i]) || isnan(dd[i])) && (isinf(cc[i]) || isnan(cc[i])))
    continue;
   assert(fabsf(dd[i] - cc[i]) < epsilon);
  }
 }
 
}

void oxb_matrix_test() {
 printf("[mat3x3] test...\n");
 mat3x3_test();
 printf("\n");
}

#endif

