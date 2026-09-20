#include "zp_physics/core2d/body2d.h"
#include "zp_physics/core2d.h"
#include "zp_physics/core2d/world2d.h"
#include <assert.h>


zp_cold zp_noinline static void init_box(zp_box2d *const zp_restrict body, const zp_create_body2d *const zp_restrict data);
zp_hot zp_inline void update_aabb(zp_body2d *const zp_restrict body);


zp_cold zp_noinline void zp_body2d_init(zp_body2d *const zp_restrict body, const void *const zp_restrict data) {
 switch((*((uint16_t*)data)) & ZP_BODY_MASK_2D) {
  case ZP_BODY_BOX_2D:
   init_box(&body->_box, (const zp_create_body2d*)data);
  break;
  default :
   assert(0);
   zp_unreachable();
  break;
 }
}



zp_hot zp_inline void update_aabb(zp_body2d *const zp_restrict body) {
 switch(body->_head._flags & ZP_BODY_MASK_2D) {
  case ZP_BODY_BOX_2D:
  {
   zp_box2d *const box = &body->_box;
   zp_vec2 rv;
   zp_complex r = body->_head._rotation;

   r = zp_abs2(r);
   rv.x = zp_fma(box->_half_size.x, r.x, box->_half_size.y * r.y);
   rv.y = zp_fma(box->_half_size.x, r.y, box->_half_size.y * r.x);
   body->_head._aabb[0] = zp_sub2(body->_head._position, rv);
   body->_head._aabb[1] = zp_add2(body->_head._position, rv);
  }
  break;
  default :
   assert(0);
   zp_unreachable();
  break;
 }
}


extern void sincosf(float x, float *s, float *c);



void zp_body2d_updatev(zp_body2d *const zp_restrict body, const void *const zp_restrict world, const float dt) {

 const zp_world2d *const ctx = (zp_world2d*)world;
 const zp_vec2 dt_vec = zp_stv2(dt);
 const zp_vec2 inv_mass = zp_stv2(body->_head._inv_mass);
 
 zp_vec2 linear_acceleration = zp_fma2(body->_head._force, inv_mass, ctx->_gravity);
 body->_head._velocity = zp_fma2(linear_acceleration, dt_vec, body->_head._velocity);

 float angular_acceleration = body->_head._torque * body->_head._inv_inertia;
 body->_head._omega = zp_fma(angular_acceleration, dt, body->_head._omega);

 body->_head._torque = 0.0f;
 body->_head._force = zp_stv2(0.0f);
}



void zp_body2d_updatep(zp_body2d *const zp_restrict body, const void *const zp_restrict world, const float dt) {
 (void)world;

 zp_vec2 dt_vec = zp_stv2(dt);

 float ld = zp_exp2(dt * body->_head._linear_damping);
 body->_head._velocity.x *= ld;
 body->_head._velocity.y *= ld;
   
 body->_head._omega *= zp_exp2(dt * body->_head._angular_damping);

 const float o_epsilon = 0.001f;
 const float v_epsilon = 0.1f;
 
 float dt_a = zp_dot2(body->_head._velocity, body->_head._velocity);
 
 uint8_t update_rotation = zp_abs(body->_head._omega) > o_epsilon;
 uint8_t update_position = dt_a > v_epsilon;
 if(update_rotation) {
  float omega = body->_head._omega * dt;
  zp_complex omega_complex;
  omega_complex.x = 1.0f - 0.5f * omega * omega;
  omega_complex.y = omega;
  body->_head._rotation = zp_cmul(omega_complex, body->_head._rotation);
  body->_head._rotation = zp_unit2(body->_head._rotation);
 }
  

 if(update_position) {
  body->_head._position = zp_fma2(body->_head._velocity, dt_vec, body->_head._position);   
 }
 if(update_rotation || update_position)
  update_aabb(body);
}



/*
 statics
*/

zp_cold zp_noinline static void init_box(zp_box2d *const zp_restrict body, const zp_create_body2d *const zp_restrict data) {
 body->_head._flags = data->_flags;
 body->_half_size = data->_half_size;
 body->_head._position = data->_position;
 body->_head._velocity = data->_velocity;
 body->_head._force = data->_force;
 body->_head._rotation = zp_crotate(data->_rotation);
 body->_head._linear_damping = zp_log2(data->_linear_damping);
 body->_head._omega = data->_omega;
 body->_head._torque = data->_torque;
 body->_head._angular_damping = zp_log2(data->_angular_damping);
 body->_head._restitution = data->_restitution;
 body->_head._friction = data->_friction;

 update_aabb((zp_body2d*)body);

 switch(data->_flags & ZP_BODY_MOVEMENT_MASK_2D) {
  case ZP_BODY_MOVEMENT_DYNAMIC_2D:
  {
   float mx = body->_half_size.x + body->_half_size.x;
   float my = body->_half_size.y + body->_half_size.y;
   float mass = mx * my * data->_density;
   body->_head._inv_mass = 1.0f / mass;
   body->_head._inv_inertia = 1.0f / (mass * (mx * mx + my * my) / 12.0f);
  }
  break;
  case ZP_BODY_MOVEMENT_KINEMATIC_2D:
   body->_head._inv_mass = 0.0f;
   body->_head._inv_inertia = 0.0f;
  break;
  case ZP_BODY_MOVEMENT_STATIC_2D:
   body->_head._inv_mass = 0.0f;
   body->_head._inv_inertia = 0.0f;
  break;
  default :
   assert(0);
   zp_unreachable();
  break;
 }
}



