#include "zp_physics/core2d/manifold2d.h"
#include "zp_physics/matrix.h"
#include <string.h>
#include <assert.h>

/* make sure it is always exactly 2 bytes */
static_assert(sizeof(zp_contact2d_id) == 2, "is not 2 bytes");


void zp_manifold2d_combine(zp_manifold2d *const zp_restrict out, const zp_manifold2d *const zp_restrict a) {
 out->_body_a = a->_body_a;
 out->_body_b = a->_body_b;
 
 zp_contact2d contacts[2];
 uint8_t contact_count = 0;
 for(uint8_t i = 0; i < a->_contact_count; i++) {
  const zp_contact2d *const contact_a = a->_contacts + i;

  int8_t j = -1;
  zp_contact2d *contact_out = NULL;
  for(uint8_t k = 0; k < out->_contact_count; k++) {
   contact_out = out->_contacts + k;
   if(contact_a->_id.value == contact_out->_id.value) {
    j = k;
    break;
   }
  }
  
  if(j > -1) {
   contacts[contact_count] = out->_contacts[j];
   contact_out = contacts + contact_count;
   contact_out->_r1 = contact_a->_r1;
   contact_out->_r2 = contact_a->_r2;
   contact_out->_normal = contact_a->_normal;
   contact_out->_depth = contact_a->_depth;
   contact_count++;
  } else {
   contacts[contact_count++] = *contact_a;
  }
 }
 out->_contact_count = contact_count;
 out->_contacts[0] = contacts[0];
 out->_contacts[1] = contacts[1];
}




zp_inline zp_const zp_contact2d_id zp_swap_contact(const zp_contact2d_id a) {
 zp_contact2d_id out;
 out.in_edge1 = a.in_edge2;
 out.in_edge2 = a.in_edge1;

 out.out_edge1 = a.out_edge2;
 out.out_edge2 = a.out_edge1;
 
 return out;
}

/***************************************
           BOX VS BOX COLLISION
****************************************/



#define OBB_EDGE_RIGHT          0
#define OBB_EDGE_LEFT           1
#define OBB_EDGE_TOP            2
#define OBB_EDGE_BOTTOM         3

/*
 lookup tables
*/
static const uint8_t BOX_SIDE_AXIS_LUT[2] = {
 1, 0
};

static const uint8_t BOX_NEGATE_AXIS[4] = {
 OBB_EDGE_LEFT, OBB_EDGE_RIGHT, OBB_EDGE_BOTTOM, OBB_EDGE_TOP
};

static const uint8_t BOX_PERP_SIDE[4][2] = {
 {OBB_EDGE_TOP, OBB_EDGE_BOTTOM},
 {OBB_EDGE_BOTTOM, OBB_EDGE_TOP},
 {OBB_EDGE_LEFT, OBB_EDGE_RIGHT},
 {OBB_EDGE_RIGHT, OBB_EDGE_LEFT},

};

/*
edges 
      1
  •--------•
  |        |
2 |        | 0
  |        |
  •--------•
      3
*/
zp_noinline static void obb_get_inc_edge(const zp_vec2 half_size, const zp_mat3x3 transform, const uint8_t index, zp_vec2 edge[2]) {
 zp_vec2 aa = zp_mul2(zp_load2(transform.v[0]), half_size);
 zp_vec2 bb = zp_mul2(zp_load2(transform.v[1]), half_size);
 zp_vec2 cc = zp_load2(transform.v[2]);

 switch(index) {
  case 0: /* right */
   edge[0] = zp_add2(zp_add2(cc, aa), bb);
   edge[1] = zp_sub2(zp_add2(cc, aa), bb);
  break;
  case 1: /* top */
   edge[0] = zp_add2(zp_add2(cc, aa), bb);
   edge[1] = zp_add2(zp_sub2(cc, aa), bb);
  break;
  case 2: /* left */
   edge[0] = zp_add2(zp_sub2(cc, aa), bb);
   edge[1] = zp_sub2(zp_sub2(cc, aa), bb);
  break;
  case 3: /* bottom */
   edge[0] = zp_sub2(zp_add2(cc, aa), bb);
   edge[1] = zp_sub2(zp_sub2(cc, aa), bb);
  break;
 }
}
/*
     •
    /  \--------> incident box
  /      \
 •        •   
  \     /  |
  | \ /    |----> clip segment
  |  •     |
  •--------•
  |        |
  |        |---> reference box
  |        |
  •--------•
     
*/
zp_noinline static uint8_t obb_clip_edge(zp_vec2 edge[2], const zp_vec2 normal, const float offset, const uint16_t edge_id, zp_contact2d_id ids[2]) {
 const float da = zp_dot2(normal, edge[0]) - offset;
 const float db = zp_dot2(normal, edge[1]) - offset;
 uint8_t poly_size = 0;
 zp_vec2 out_edge[2];
 zp_contact2d_id out_ids[2];
 
 const float c_epsilon = 1e-3f;

 if(da < c_epsilon) {
 	out_edge[poly_size] = edge[0];
 	out_ids[poly_size] = ids[0];
 	poly_size++;
 }
 if(db < c_epsilon) {
  out_edge[poly_size] = edge[1];
  out_ids[poly_size] = ids[1];
  poly_size++;
 }
 
 if((da * db) < 0.0f && poly_size < 2) {
  zp_vec2 t = zp_stv2(db / (db - da));
  out_edge[poly_size] = zp_fma2(zp_sub2(edge[0], edge[1]), t, edge[1]);

  zp_contact2d_id intersection_id;
  if(da > 0.0f) {
   intersection_id.value = ids[0].value;
   intersection_id.in_edge1 = edge_id;
   intersection_id.in_edge2 = 0xF;
  } else {
   intersection_id.value = ids[1].value;
   intersection_id.out_edge1 = edge_id;
   intersection_id.out_edge2 = 0xF;
  }
  out_ids[poly_size] = intersection_id;
		poly_size++;
 }
 assert(poly_size <= 2); 
 edge[0] = out_edge[0];
 edge[1] = out_edge[1];
 ids[0] = out_ids[0];
 ids[1] = out_ids[1];
 return poly_size;
}

/*
 Real-Time Collision Detection by Christer Ericson
 page 143
 
 OBB vs OBB
 src : http://www.r-5.org/files/books/computers/algo-list/realtime-3d/Christer_Ericson-Real-Time_Collision_Detection-EN.pdf
*/
zp_inline float obb_extent_size(const zp_mat3x3 obj_t, const zp_vec2 obj_hs, const zp_vec2 axis) {
 return zp_fma(obj_hs.x, zp_abs(zp_dot2(axis, zp_load2(obj_t.v[0]))), obj_hs.y * zp_abs(zp_dot2(axis, zp_load2(obj_t.v[1]))) );
}

/*
 instead of polygon, box treated as obb.
*/
uint8_t zp_manifold2d_box_vs_box(zp_manifold2d *const zp_restrict out, const zp_box2d *const zp_restrict a, const zp_box2d *const zp_restrict b) {
 zp_vec2 half_size[2];
 zp_complex rotation_a = a->_head._rotation;
 zp_vec2 position_a = a->_head._position;
 half_size[0] = a->_half_size;
 
 zp_compiler_memory_barrier();

 zp_complex rotation_b = b->_head._rotation;
 zp_vec2 position_b = b->_head._position;
 half_size[1] = b->_half_size;
 
 zp_compiler_memory_barrier();
 
 zp_mat3x3 transform_a = zp_mat3x3_transform(rotation_a, position_a);
 zp_mat3x3 transform_b = zp_mat3x3_transform(rotation_b, position_b);
 
 zp_vec2 diff_pos = zp_sub2(position_b, position_a);
 
 float ra, rb, dist, depth_a, depth_b;
 uint8_t index_a, index_b;
 zp_vec2 axis;
 
 depth_a = zp_inf();
 depth_b = zp_inf();
 
 /*
  fast SAT collision without using vertices.
 */
 for(uint8_t i = 0; i < 2; i++) {
  axis = zp_load2(transform_a.v[i]);
  dist = zp_abs(zp_dot2(diff_pos, axis));   
  ra = half_size[0].arr[i];
  rb = obb_extent_size(transform_b, half_size[1], axis);
  if(dist > (ra + rb)) return 0;
  dist = ra + rb - dist;
  if(dist < depth_a) {
   depth_a = dist;
   index_a = i;
  }
 }

 for(uint8_t i = 0; i < 2; i++) {
  axis = zp_load2(transform_b.v[i]);
  dist = zp_abs(zp_dot2(diff_pos, axis));   
  ra = obb_extent_size(transform_a, half_size[0], axis);
  rb = half_size[1].arr[i];
  if(dist > (ra + rb)) return 0;
  dist = ra + rb - dist;
  if(dist < depth_b) {
   depth_b = dist;
   index_b = i;
  }
 }
 

 uint8_t ref_index, inc_index;
 zp_mat3x3 ref_transform, inc_transform;
 zp_vec2 ref_half_size, inc_half_size;
 uint8_t swapped;
 const zp_box2d *ref_body, *inc_body;

 swapped = depth_b < depth_a;

 zp_contact2d_id ids[2];
 ids[0].value = 0;
 ids[1].value = 0;
 if(swapped) {
  ref_index = index_b;
  inc_index = index_a;
  ref_transform = transform_b;
  inc_transform = transform_a;
  ref_half_size = half_size[1];
  inc_half_size = half_size[0];
  ref_body = b;
  inc_body = a;
  
  diff_pos = zp_neg2(diff_pos);
 } else {
  ref_index = index_a;
  inc_index = index_b;
  ref_transform = transform_a;
  inc_transform = transform_b;
  ref_half_size = half_size[0];
  inc_half_size = half_size[1];
  ref_body = a;
  inc_body = b;
 }


 zp_vec2 ref_axis = zp_load2(ref_transform.v[ref_index]);

 uint16_t tmp_ref_index = (uint16_t)ref_index;

 if(zp_dot2(diff_pos, ref_axis) < 0.0f) {
  ref_axis = zp_neg2(ref_axis);
  tmp_ref_index += 2;
 }


 ids[0].in_edge1 = BOX_PERP_SIDE[tmp_ref_index][0];
 ids[0].out_edge1 = tmp_ref_index;

 ids[1].in_edge1 = tmp_ref_index;
 ids[1].out_edge1 = BOX_PERP_SIDE[tmp_ref_index][1];



 uint8_t flip = zp_dot2(ref_axis, zp_load2(inc_transform.v[inc_index])) > 0.0f;
 uint8_t inc_edge_index = inc_index + (flip ? 2 : 0);
 
 ids[0].in_edge2 = BOX_PERP_SIDE[inc_edge_index][0];
 ids[0].out_edge2 = inc_edge_index;

 ids[1].in_edge2 = inc_edge_index;
 ids[1].out_edge2 = BOX_PERP_SIDE[inc_edge_index][1];

 zp_vec2 edge[2];
 obb_get_inc_edge(inc_half_size, inc_transform, inc_edge_index, edge);

 uint8_t vertex_count;
 const uint8_t ref_side_index = BOX_SIDE_AXIS_LUT[ref_index];
 zp_vec2 ref_side_axis = zp_load2(ref_transform.v[ref_side_index]);
 float pos_offset = zp_dot2(zp_load2(ref_transform.v[2]), ref_side_axis) + ref_half_size.arr[ref_side_index];
 float neg_offset = -zp_dot2(zp_load2(ref_transform.v[2]), ref_side_axis) + ref_half_size.arr[ref_side_index];

 vertex_count = obb_clip_edge(edge, zp_neg2(ref_side_axis), neg_offset, BOX_NEGATE_AXIS[ref_side_index], ids);
 if(vertex_count < 2) return 0;

 vertex_count = obb_clip_edge(edge, ref_side_axis, pos_offset, ref_side_index, ids); 
 if(vertex_count < 2) return 0;

 
 const float front = zp_dot2(zp_load2(ref_transform.v[2]), ref_axis) + ref_half_size.arr[ref_index];
 if(swapped) {

  /* expected to be normalized, so the inverse is just as simple as conjugate */ 
  zp_complex inv_ref_rotation = zp_cconj(rotation_b);
  zp_complex inv_inc_rotation = zp_cconj(rotation_a);

  for(uint8_t i = 0; i < 2; i++) {
   float depth = zp_dot2(ref_axis, edge[i]) - front;
   if(depth < 0.0f) {
    zp_contact2d *const contact = out->_contacts + out->_contact_count;
    contact->_id = zp_swap_contact(ids[i]);
    zp_vec2 local_ra = zp_sub2(edge[i], zp_load2(ref_transform.basis[2].arr));
    zp_vec2 local_rb = zp_sub2(edge[i], zp_load2(inc_transform.basis[2].arr));
    contact->_r1 = zp_cmul(inv_ref_rotation, zp_as_complex2(local_ra));
    contact->_r2 = zp_cmul(inv_inc_rotation, zp_as_complex2(local_rb));
    contact->_depth = depth;
    contact->_normal = ref_axis;
    out->_contact_count++;
   }
  }
 } else {

  /* expected to be normalized, so the inverse is just as simple as conjugate */ 
  zp_complex inv_ref_rotation = zp_cconj(rotation_a);
  zp_complex inv_inc_rotation = zp_cconj(rotation_b);

  for(uint8_t i = 0; i < 2; i++) {
   float depth = zp_dot2(ref_axis, edge[i]) - front;
   if(depth < 0.0f) {
    zp_contact2d *const contact = out->_contacts + out->_contact_count;
    contact->_id = ids[i];
    zp_vec2 local_ra = zp_sub2(edge[i], zp_load2(ref_transform.basis[2].arr));
    zp_vec2 local_rb = zp_sub2(edge[i], zp_load2(inc_transform.basis[2].arr));
    contact->_r1 = zp_cmul(inv_ref_rotation, zp_as_complex2(local_ra));
    contact->_r2 = zp_cmul(inv_inc_rotation, zp_as_complex2(local_rb));
    contact->_depth = depth;
    contact->_normal = ref_axis;
    out->_contact_count++;
   }
  }
 }
 out->_body_a = (zp_container_id)ref_body->_head._id;
 out->_body_b = (zp_container_id)inc_body->_head._id;
 return out->_contact_count;
}



