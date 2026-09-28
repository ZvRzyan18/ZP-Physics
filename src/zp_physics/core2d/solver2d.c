#include "zp_physics/core2d/solver2d.h"
#include "zp_physics/core2d/world2d.h"
#include "zp_physics/core2d.h"

/*
 src : https://github.com/erincatto/solver2d
 Temporal Gauss Seidel + soft constraints
*/


void zp_manifold2d_soft_prepare_contact(zp_manifold2d *const zp_restrict m, void *const zp_restrict w, const zp_solver_input2d *const input) {
 zp_world2d *const world = (zp_world2d*)w;

 zp_body2d *const body_a = (zp_body2d*)zp_container_get(&world->_body_container, m->_body_a);
 zp_body2d *const body_b = (zp_body2d*)zp_container_get(&world->_body_container, m->_body_b);
/*
 if(zp_body2d_asleep(body_a) && zp_body2d_asleep(body_b))
  return;
  */
 zp_compiler_memory_barrier();

 zp_vec2 a_pos = body_a->_head._position;
 zp_complex a_rot = body_a->_head._rotation;
 zp_vec2 a_velocity = body_a->_head._velocity;
 float a_inv_mass = body_a->_head._inv_mass;
 float a_inv_inertia = body_a->_head._inv_inertia;
 float a_restitution = body_a->_head._restitution;
 float a_omega = body_a->_head._omega;

 zp_compiler_memory_barrier();

 zp_vec2 b_pos = body_b->_head._position;
 zp_complex b_rot = body_b->_head._rotation;
 zp_vec2 b_velocity = body_b->_head._velocity;
 float b_inv_mass = body_b->_head._inv_mass;
 float b_inv_inertia = body_b->_head._inv_inertia;
 float b_restitution = body_b->_head._restitution;
 float b_omega = body_b->_head._omega;



 (void)input;

 float inv_mass = a_inv_mass + b_inv_mass;
 float e = zp_min(a_restitution, b_restitution);

 
 for(uint8_t i = 0; i < m->_contact_count; i++) {
 	zp_contact2d *const contact = m->_contacts + i;

  contact->_updated_r1 = zp_as_vector2(zp_cmul(a_rot, contact->_r1));
  contact->_updated_r2 = zp_as_vector2(zp_cmul(b_rot, contact->_r2));

  zp_vec2 normal = contact->_normal;
	 zp_vec2 tangent = zp_perp2(normal);

  zp_vec2 r1 = contact->_updated_r1;
  zp_vec2 r2 = contact->_updated_r2;
  
  zp_vec2 pa = zp_add2(a_pos, r1);
  zp_vec2 pb = zp_add2(b_pos, r2);

  contact->_target_depth = zp_dot2(zp_sub2(pb, pa), normal) + contact->_depth;
  
  float r1_length_squared = zp_dot2(r1, r1);
  float r2_length_squared = zp_dot2(r2, r2);
 
  {
	  float rn1_squared = zp_dot2(r1, normal);
	  float rn2_squared = zp_dot2(r2, normal);
	  rn1_squared *= rn1_squared;
	  rn2_squared *= rn2_squared;
	 
		 float k_normal = zp_fma(zp_fma(a_inv_inertia, (r1_length_squared - rn1_squared), b_inv_inertia), (r2_length_squared - rn2_squared), inv_mass);
  	contact->_mass_normal = 1.0f / k_normal;
  }

  {
	  float rt1_squared = zp_dot2(r1, tangent);
	  float rt2_squared = zp_dot2(r2, tangent);
	  rt1_squared *= rt1_squared;
	  rt2_squared *= rt2_squared;
	 
	  float k_tangent = zp_fma(zp_fma(a_inv_inertia, (r1_length_squared - rt1_squared), b_inv_inertia), (r2_length_squared - rt2_squared), inv_mass);
	  contact->_mass_tangent = 1.0f / k_tangent;
  }
  /*
   Not correct. but by allowing only the linear relative velocity
   makes it more predictable.
  */
  zp_vec2 relative_vel = zp_sub2(b_velocity, a_velocity);
  float impact_speed = zp_dot2(relative_vel, normal);
  /* 
   hyperbolic tangent gives a smooth transition betweew low and max resitution coeffs


   tanh(x)
   1.0 ---------------------------------
                      ____
                  ___
              ___
           __
         __
       _
     -
    -
  -------------------------------------

  conditional by threshold
  -------------------------------------
                  |-------------------
                  |
                  |
                  |
                  |
                  |
   _______________
  -------------------------------------
  
  */
  float x = impact_speed * 0.2f;
  contact->_bias = (impact_speed * e) * zp_tanh(zp_min(x, 0.0f));
 }

 float la = zp_dot2(a_velocity, a_velocity);
 float lb = zp_dot2(b_velocity, b_velocity);
 float oa = zp_abs(a_omega);
 float ob = zp_abs(b_omega);

 {
 const float vt =  2.2f;
 const float ot =  0.8f;
 if(la > vt || oa > ot)
  zp_body2d_awake(body_a);
 if(lb > vt || ob > ot)
  zp_body2d_awake(body_b);
 }
}



void zp_manifold2d_soft_presolve_contact(zp_manifold2d *const zp_restrict m, void *const zp_restrict w, const zp_solver_input2d *const input) {
 (void)input;
 
 zp_world2d *const world = (zp_world2d*)w;

 zp_body2d *const body_a = (zp_body2d*)zp_container_get(&world->_body_container, m->_body_a);
 zp_body2d *const body_b = (zp_body2d*)zp_container_get(&world->_body_container, m->_body_b);
/*
 if(zp_body2d_asleep(body_a) && zp_body2d_asleep(body_b))
  return;
*/
 zp_compiler_memory_barrier();

 zp_vec2 a_velocity = body_a->_head._velocity;
 float a_omega = body_a->_head._omega;
 float a_inv_mass = body_a->_head._inv_mass;
 float a_inv_inertia = body_a->_head._inv_inertia;

 zp_compiler_memory_barrier();

 zp_vec2 b_velocity = body_b->_head._velocity;
 float b_omega = body_b->_head._omega;
 float b_inv_mass = body_b->_head._inv_mass;
 float b_inv_inertia = body_b->_head._inv_inertia;

 for(uint8_t i = 0; i < m->_contact_count; i++) {
 	zp_contact2d *const contact = m->_contacts + i;

  zp_vec2 r1 = contact->_updated_r1;
  zp_vec2 r2 = contact->_updated_r2;

  zp_vec2 normal = contact->_normal;
	 zp_vec2 tangent = zp_perp2(normal);

	 zp_vec2 impulse = zp_add2(zp_mul2(zp_stv2(contact->_accumulated_normal), normal), zp_mul2(zp_stv2(contact->_accumulated_tangent), tangent));
 	a_velocity = zp_sub2(a_velocity, zp_mul2(zp_stv2(a_inv_mass), impulse));
	 a_omega -= a_inv_inertia * zp_cross2(r1, impulse);
		 
	 b_velocity = zp_add2(b_velocity, zp_mul2(zp_stv2(b_inv_mass), impulse));
	 b_omega += b_inv_inertia * zp_cross2(r2, impulse);
 } 
 body_a->_head._velocity = a_velocity;
 body_a->_head._omega = a_omega;

 zp_compiler_memory_barrier();

 body_b->_head._velocity = b_velocity;
 body_b->_head._omega = b_omega;
}


void zp_manifold2d_soft_solve_contact(zp_manifold2d *const zp_restrict m, void *const zp_restrict w, const zp_solver_input2d *const input) {
 zp_world2d *const world = (zp_world2d*)w;
 
 zp_body2d *const body_a = (zp_body2d*)zp_container_get(&world->_body_container, m->_body_a);
 zp_body2d *const body_b = (zp_body2d*)zp_container_get(&world->_body_container, m->_body_b);
/*
 if(zp_body2d_asleep(body_a) && zp_body2d_asleep(body_b))
  return;
*/
 zp_compiler_memory_barrier();

 zp_vec2 a_velocity = body_a->_head._velocity;
 float a_omega = body_a->_head._omega;
 float a_inv_mass = body_a->_head._inv_mass;
 float a_inv_inertia = body_a->_head._inv_inertia;
 float a_friction = body_a->_head._friction;

 zp_compiler_memory_barrier();

 zp_vec2 b_velocity = body_b->_head._velocity;
 float b_omega = body_b->_head._omega;
 float b_inv_mass = body_b->_head._inv_mass;
 float b_inv_inertia = body_b->_head._inv_inertia;
 float b_friction = body_b->_head._friction;


 float friction = zp_sqrt(a_friction * b_friction); 
 
 zp_vec2 r1, r2;
 zp_vec2 va, vb;
 zp_vec2 relative_vel;
 zp_vec2 impulse;
 float j, depth;
 

 for(uint8_t i = 0; i < m->_contact_count; i++) {
 	zp_contact2d *const contact = m->_contacts + i;
  
  r1 = contact->_updated_r1;
  r2 = contact->_updated_r2;
  depth = contact->_target_depth;

  zp_vec2 normal = contact->_normal;

	 va = zp_add2(a_velocity, zp_cross_sv2(a_omega, r1));
  vb = zp_add2(b_velocity, zp_cross_sv2(b_omega, r2));
	 relative_vel = zp_sub2(vb, va);

		float vn = zp_dot2(relative_vel, normal);


  float slop = 0.004f;
  float penetration_error = 0.0f;
  
  if(depth < 0.0f) {
   float error = zp_min((depth + slop), 0.0f);
   penetration_error = -error * input->_bias_ratio;
  } else {
   penetration_error = depth * input->_inv_dt;  
  }

  /* logarithmic with damping-like behaviour */
  penetration_error = zp_log2(penetration_error + 1.0f) * 0.69314718f;
  penetration_error = -(vn - contact->_bias) + penetration_error;
  
  {
   float mass = contact->_mass_normal * input->_mass_coeff;
	  j = mass * penetration_error;
   j -= input->_impulse_coeff * contact->_accumulated_normal;
  }
  {
		 float accumulated_normal = contact->_accumulated_normal;
	 	contact->_accumulated_normal = zp_max(accumulated_normal + j, 0.0f);
		 j = contact->_accumulated_normal - accumulated_normal;
  }
  
	 impulse = zp_mul2(zp_stv2(j), normal);
	 
 	a_velocity = zp_sub2(a_velocity, zp_mul2(zp_stv2(a_inv_mass), impulse));
	 a_omega -= a_inv_inertia * zp_cross2(r1, impulse);
		 
	 b_velocity = zp_add2(b_velocity, zp_mul2(zp_stv2(b_inv_mass), impulse));
	 b_omega += b_inv_inertia * zp_cross2(r2, impulse);


 /*
  tangent / friction impulse
  
   ^
  /|\
   |
   | --- normal direction
   | 
   |       ______ tangent direction
   |       |
   •------------- >
  
  • resist sliding behaviour 
  • always perpendicular or 90°
 */
  zp_vec2 tangent = zp_perp2(normal);
  
  va = zp_add2(a_velocity, zp_cross_sv2(a_omega, r1));
  vb = zp_add2(b_velocity, zp_cross_sv2(b_omega, r2));
	 
	 relative_vel = zp_sub2(vb, va);
	 
		float vt = zp_dot2(relative_vel, tangent);
 	j = contact->_mass_tangent * -vt;


		float friction_range = friction * contact->_accumulated_normal;

  {
	 	float prev_tangent = contact->_accumulated_tangent;
	 	contact->_accumulated_tangent = zp_min(zp_max(prev_tangent + j, -friction_range), friction_range);
		 j = contact->_accumulated_tangent - prev_tangent;
  }
  
  impulse = zp_mul2(zp_stv2(j), tangent);
 	a_velocity = zp_sub2(a_velocity, zp_mul2(zp_stv2(a_inv_mass), impulse));
	 a_omega -= a_inv_inertia * zp_cross2(r1, impulse);
		 
	 b_velocity = zp_add2(b_velocity, zp_mul2(zp_stv2(b_inv_mass), impulse));
	 b_omega += b_inv_inertia * zp_cross2(r2, impulse);
 }
 body_a->_head._velocity = a_velocity;
 body_a->_head._omega = a_omega;

 zp_compiler_memory_barrier();

 body_b->_head._velocity = b_velocity;
 body_b->_head._omega = b_omega;
}







