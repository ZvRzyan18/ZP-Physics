#include "oxb1tz/core2d/manifold2d.h"
#include "oxb1tz/matrix.h"

static_assert(sizeof(oxb_contact2d_id) == 2, "not 2 bytes");


void oxb_manifold2d_combine(oxb_manifold2d *oxb_restrict const out, const oxb_manifold2d *oxb_restrict const m) {
 out->_a = m->_a;
 out->_b = m->_b;
 
 oxb_contact2d contacts[2];
 uint8_t contact_count = 0;
 for(size_t i = 0; i < m->_contact_count; i++) {
  const oxb_contact2d *const contact_a = m->_contacts + i;

  int32_t j = -1;
  oxb_contact2d *contact_out = NULL;
  for(size_t k = 0; k < out->_contact_count; k++) {
   contact_out = out->_contacts + k;
   if(contact_a->_id._value == contact_out->_id._value) {
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



oxb_inline oxb_const oxb_contact2d_id swap_contact(const oxb_contact2d_id a) {
 oxb_contact2d_id out;
 out._in_edge1 = a._in_edge2;
 out._in_edge2 = a._in_edge1;

 out._out_edge1 = a._out_edge2;
 out._out_edge2 = a._out_edge1;
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


static void obb_get_inc_edge(const oxb_vec2 half_size, const oxb_mat3x3 transform, const uint8_t index, oxb_vec2 edge[2]) {
 oxb_vec2 scale_x = oxb_vec2_splat(oxb_vec2_getx(half_size));
 oxb_vec2 scale_y = oxb_vec2_splat(oxb_vec2_gety(half_size));
 
 oxb_vec2 aa = oxb_vec3_low(transform._basis[0]);
 oxb_vec2 bb = oxb_vec3_low(transform._basis[1]);
 oxb_vec2 cc = oxb_vec3_low(transform._basis[2]);
 
 aa = oxb_vec2_mul(scale_x, aa);
 bb = oxb_vec2_mul(scale_y, bb);
 
 switch(index) {
  case 0: /* right */
   edge[0] = oxb_vec2_add(oxb_vec2_add(cc, aa), bb); /* top-right */ 
   edge[1] = oxb_vec2_sub(oxb_vec2_add(cc, aa), bb); /* bottom-right */
  break;
  case 1: /* top */
   edge[0] = oxb_vec2_add(oxb_vec2_sub(cc, aa), bb); /* top-left */
   edge[1] = oxb_vec2_add(oxb_vec2_add(cc, aa), bb); /* top-right */
  break;
  case 2: /* left */
   edge[0] = oxb_vec2_sub(oxb_vec2_sub(cc, aa), bb); /* bottom-left */
   edge[1] = oxb_vec2_add(oxb_vec2_sub(cc, aa), bb); /* top-left */
  break;
  case 3: /* bottom */
   edge[0] = oxb_vec2_sub(oxb_vec2_add(cc, aa), bb); /* bottom-right */
   edge[1] = oxb_vec2_sub(oxb_vec2_sub(cc, aa), bb); /* bottom-left */
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
static uint8_t obb_clip_edge(oxb_vec2 edge[2], const oxb_vec2 normal, const float offset, const uint8_t edge_id, oxb_contact2d_id ids[2]) {
 const float da = oxb_vec2_dot(normal, edge[0]) - offset;
 const float db = oxb_vec2_dot(normal, edge[1]) - offset;
 uint8_t poly_size = 0;
 oxb_vec2 out_edge[2];
 oxb_contact2d_id out_ids[2];
 
 const float c_epsilon = 1e-3f;

 if (da < c_epsilon) {
  out_edge[poly_size] = edge[0];
  out_ids[poly_size] = ids[0];
  poly_size++;
 }
 if (db < c_epsilon) {
  out_edge[poly_size] = edge[1];
  out_ids[poly_size] = ids[1];
  poly_size++;
 }
 
 if ((da * db) < 0.0f && poly_size < 2) {
  oxb_vec2 t = oxb_vec2_splat(db / (db - da));
  out_edge[poly_size] = oxb_vec2_fma(oxb_vec2_sub(edge[0], edge[1]), t, edge[1]);

  oxb_contact2d_id intersection_id;
  if (da > 0.0f) {
   intersection_id = ids[0];
   intersection_id._in_edge1 = edge_id;
   intersection_id._in_edge2 = ids[0]._in_edge2; 
  } else {
   intersection_id = ids[1];
   intersection_id._out_edge1 = edge_id;
   intersection_id._out_edge2 = ids[1]._out_edge2;
  }
  out_ids[poly_size] = intersection_id;
  poly_size++;
 }
 
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
oxb_inline float obb_extent_size(const oxb_mat3x3 obj_t, const oxb_vec2 obj_hs, const oxb_vec2 axis) {
 float a = oxb_vec2_dot(axis, oxb_vec3_low(obj_t._basis[0]));
 float b = oxb_vec2_dot(axis, oxb_vec3_low(obj_t._basis[1]));
 oxb_vec2 tmp = oxb_vec2_make(a, b);
 return oxb_vec2_dot(obj_hs, oxb_vec2_abs(tmp));
}


/*
 instead of polygon, box treated as obb.
*/
bool oxb_manifold2d_box_vs_box(oxb_manifold2d *oxb_restrict const out, const oxb_collider2d *oxb_restrict const a, const oxb_collider2d *oxb_restrict const b) {
 /* aabb early rejection */
 if((a->_aabb._max._arr[0] < b->_aabb._min._arr[0] || a->_aabb._min._arr[0] > b->_aabb._max._arr[0] || a->_aabb._max._arr[1] < b->_aabb._min._arr[1] || a->_aabb._min._arr[1] > b->_aabb._max._arr[1]))
  return false;
 
 oxb_vec2 half_size[2];
 
 oxb_complex rotation_a = oxb_vec2_load(a->_rotation._arr);
 oxb_vec2 position_a = oxb_vec2_load(a->_position._arr);
 half_size[0] = oxb_vec2_load(a->_shape._box._half_size._arr);
 
 oxb_complex rotation_b = oxb_vec2_load(b->_rotation._arr);
 oxb_vec2 position_b = oxb_vec2_load(b->_position._arr);
 half_size[1] = oxb_vec2_load(b->_shape._box._half_size._arr);
 
   
 oxb_mat3x3 transform_a = oxb_mat3x3_transform(rotation_a, position_a);
 oxb_mat3x3 transform_b = oxb_mat3x3_transform(rotation_b, position_b);
 
 oxb_vec2 diff_pos = oxb_vec2_sub(position_b, position_a);
 
 float ra, rb, dist, depth_a, depth_b;
 uint8_t index_a, index_b;
 oxb_vec2 axis;
 
 depth_a = oxb_inf();
 depth_b = oxb_inf();

 
 /*
  fast SAT collision without using vertices.
  hope that compiler unroll this.
 */
 for(size_t i = 0; i < 2; i++) {
  axis = oxb_vec3_low(transform_a._basis[i]);
  dist = oxb_abs(oxb_vec2_dot(diff_pos, axis));   
  ra = oxb_vec2_get(half_size[0], i);
  rb = obb_extent_size(transform_b, half_size[1], axis);
  if(dist > (ra + rb)) return false;
  dist = ra + rb - dist;
  if(dist < depth_a) {
   depth_a = dist;
   index_a = i;
  }
 }

 for(size_t i = 0; i < 2; i++) {
  axis = oxb_vec3_low(transform_b._basis[i]);
  dist = oxb_abs(oxb_vec2_dot(diff_pos, axis));   
  ra = obb_extent_size(transform_a, half_size[0], axis);
  rb = oxb_vec2_get(half_size[1], i);
  if(dist > (ra + rb)) return false;
  dist = ra + rb - dist;
  if(dist < depth_b) {
   depth_b = dist;
   index_b = i;
  }
 }
 
 /* sensor collision, no need manifold creation. */
 if(a->_allocation._type == 3 || b->_allocation._type == 3)
  return true;
 
 /* manifold generation */
 uint8_t ref_index, inc_index;
 oxb_mat3x3 ref_transform, inc_transform;
 oxb_vec2 ref_half_size, inc_half_size;
 uint8_t swapped;
 const oxb_collider2d *ref_collider;
 const oxb_collider2d *inc_collider;

 swapped = depth_b < depth_a;

 oxb_contact2d_id ids[2];
 ids[0]._value = 0;
 ids[1]._value = 0;
 if(swapped) {
  ref_index = index_b;
  inc_index = index_a;
  ref_transform = transform_b;
  inc_transform = transform_a;
  ref_half_size = half_size[1];
  inc_half_size = half_size[0];
  ref_collider = b;
  inc_collider = a;
  
  diff_pos = oxb_vec2_neg(diff_pos);
 } else {
  ref_index = index_a;
  inc_index = index_b;
  ref_transform = transform_a;
  inc_transform = transform_b;
  ref_half_size = half_size[0];
  inc_half_size = half_size[1];
  ref_collider = a;
  inc_collider = b;
 }



 oxb_vec2 ref_axis = oxb_vec3_low(ref_transform._basis[ref_index]);
 uint16_t tmp_ref_index = (uint16_t)ref_index;

 if(oxb_vec2_dot(diff_pos, ref_axis) < 0.0f) {
  ref_axis = oxb_vec2_neg(ref_axis);
  tmp_ref_index += 2;
 }


 ids[0]._in_edge1 = BOX_PERP_SIDE[tmp_ref_index][0];
 ids[0]._out_edge1 = tmp_ref_index;

 ids[1]._in_edge1 = tmp_ref_index;
 ids[1]._out_edge1 = BOX_PERP_SIDE[tmp_ref_index][1];


 float d0 = oxb_vec2_dot(ref_axis, oxb_vec3_low(inc_transform._basis[0]));
 float d1 = oxb_vec2_dot(ref_axis, oxb_vec3_low(inc_transform._basis[1]));
 uint8_t flip;
 if(oxb_abs(d0) > oxb_abs(d1)) {
  inc_index = 0;
  flip = d0 > 0.0f;
 } else {
  inc_index = 1;
  flip = d1 > 0.0f;
 }
 uint8_t inc_edge_index = inc_index + (flip ? 2 : 0);

 
 ids[0]._in_edge2 = BOX_PERP_SIDE[inc_edge_index][0];
 ids[0]._out_edge2 = inc_edge_index;

 ids[1]._in_edge2 = inc_edge_index;
 ids[1]._out_edge2 = BOX_PERP_SIDE[inc_edge_index][1];

 oxb_vec2 edge[2];
 obb_get_inc_edge(inc_half_size, inc_transform, inc_edge_index, edge);


 uint8_t vertex_count;
 const uint8_t ref_side_index = BOX_SIDE_AXIS_LUT[ref_index];
 oxb_vec2 ref_side_axis = oxb_vec3_low(ref_transform._basis[ref_side_index]);
 float pos_offset = oxb_vec2_dot(oxb_vec3_low(ref_transform._basis[2]), ref_side_axis) + oxb_vec2_get(ref_half_size, ref_side_index);
 float neg_offset = -oxb_vec2_dot(oxb_vec3_low(ref_transform._basis[2]), ref_side_axis) + oxb_vec2_get(ref_half_size, ref_side_index);

 vertex_count = obb_clip_edge(edge, oxb_vec2_neg(ref_side_axis), neg_offset, BOX_NEGATE_AXIS[ref_side_index], ids);
 vertex_count = obb_clip_edge(edge, ref_side_axis, pos_offset, ref_side_index, ids); 


 const float front = oxb_vec2_dot(oxb_vec3_low(ref_transform._basis[2]), ref_axis) + oxb_vec2_get(ref_half_size, ref_index);
 if(swapped) {

  for(uint8_t i = 0; i < 2; i++) {
   float depth = oxb_vec2_dot(ref_axis, edge[i]) - front;
   if(depth < 0.0f) {
    oxb_contact2d *const contact = out->_contacts + out->_contact_count;
    contact->_id = swap_contact(ids[i]);
    oxb_vec2_store(contact->_r1._arr, oxb_vec2_sub(edge[i], oxb_vec3_low(ref_transform._basis[2])));
    oxb_vec2_store(contact->_r2._arr, oxb_vec2_sub(edge[i], oxb_vec3_low(inc_transform._basis[2])));
    contact->_depth = depth;
    oxb_vec2_store(contact->_normal._arr, ref_axis);
    out->_contact_count++;
   }
  }
  
 } else {

  for(uint8_t i = 0; i < 2; i++) {
   float depth = oxb_vec2_dot(ref_axis, edge[i]) - front;
   if(depth < 0.0f) {
    oxb_contact2d *const contact = out->_contacts + out->_contact_count;
    contact->_id = ids[i];
    oxb_vec2_store(contact->_r1._arr, oxb_vec2_sub(edge[i], oxb_vec3_low(ref_transform._basis[2])));
    oxb_vec2_store(contact->_r2._arr, oxb_vec2_sub(edge[i], oxb_vec3_low(inc_transform._basis[2])));
    contact->_depth = depth;
    oxb_vec2_store(contact->_normal._arr, ref_axis);
    out->_contact_count++;
   }
  }
  
 }
 out->_a = ref_collider->_allocation;
 out->_b = inc_collider->_allocation;
 
 return out->_contact_count > 0;
}

