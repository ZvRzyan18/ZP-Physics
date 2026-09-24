#ifndef ZP_WORLD2D_H
#define ZP_WORLD2D_H

#include "zp_physics/core2d.h"
#include "zp_physics/container.h"
#include "zp_physics/core2d/contacthash2d.h"
#include "zp_physics/core2d/broadphase2d.h"


struct zp_world2d {
 zp_container     _body_container;
 zp_contacthash2d _contacts;
 zp_broadphase2d  _broadphase;
 
 zp_vec2          _gravity;

 size_t           _dynamic_count;
 size_t           _kinematic_count;

 float            _inv_timestep_substep;
 float            _growth_base;
 
 float            _contact_constraint_hertz;
 float            _contact_constraint_damping_ratio;
 float            _contact_constraint_omega;

 uint16_t         _solver_substeps;
 uint16_t         _time_substeps;
  
 uint8_t          _sort_body;
 uint8_t          _broadphase_frame_tolerance;

 uint8_t          _dt_limit;
};

#endif

