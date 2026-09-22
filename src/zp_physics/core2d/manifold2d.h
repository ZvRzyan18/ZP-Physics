#ifndef ZP_MANIFOLD2D_H
#define ZP_MANIFOLD2D_H

#include "zp_physics/vector.h"
#include "zp_physics/core2d/body2d.h"

/*
 feature pair of contact id.
 this uses to be able to make contact persistent across multiple frames
 and make it able to reorder the contacts as needed.
 
 it uses bitfields to minimize memory consumption,
 since 2 ^ 4 = 16, its more than enough for the ids
*/
typedef union {
 struct {
  unsigned char in_edge1  : 4;
  unsigned char out_edge1 : 4;
  unsigned char in_edge2  : 4;
  unsigned char out_edge2 : 4;
 };
 uint16_t value;
} zp_contact2d_id;
/*
 per contact points info
*/
typedef struct {
 zp_complex  _r1;
 zp_complex  _r2;
 
 zp_vec2  _updated_r1;
 zp_vec2  _updated_r2;

 zp_vec2  _normal;
 float    _depth;
 float    _updated_depth;
 
 float    _mass_normal;
 float    _mass_tangent;
 float    _bias;
  
 float    _accumulated_normal;
 float    _accumulated_tangent;
 zp_contact2d_id _id;
} zp_contact2d;

typedef struct {
 zp_contact2d    _contacts[2];
 zp_container_id _body_a;
 zp_container_id _body_b;
 
 uint8_t         _contact_count;
} zp_manifold2d;



void zp_manifold2d_combine(zp_manifold2d *const zp_restrict out, const zp_manifold2d *const zp_restrict a);
uint8_t zp_manifold2d_box_vs_box(zp_manifold2d *const zp_restrict out, const zp_box2d *const zp_restrict a, const zp_box2d *const zp_restrict b);


#endif

