#ifndef ZP_WORLD2D_H
#define ZP_WORLD2D_H

#include "zp_physics/core2d.h"
#include "zp_physics/container.h"
#include "zp_physics/core2d/contacthash2d.h"

#include "zp_physics/core2d/broadphase2d.h"


struct zp_world2d {
 /* NOTE : this can change anytime, breaking compatibility, so do not access */
 zp_container     _body_container;
 zp_contacthash2d _contacts;
 zp_broadphase2d  _broadphase;
 
 zp_vec2          _gravity;
 float            _inv_timestep_substep;
 float            _growth_base;
 uint16_t         _solver_substeps;
 uint16_t         _time_substeps;
 
 zp_container_id  _dynamic_count;
 zp_container_id  _kinematic_count;
 uint8_t          _sort_body;
};

#endif

