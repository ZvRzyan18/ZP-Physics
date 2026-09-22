#ifndef ZP_BODY2D_H
#define ZP_BODY2D_H

#include "zp_physics/core2d.h"
#include "zp_physics/complex.h"
#include "zp_physics/container.h"
#include "zp_physics/pool.h"


typedef struct {
 zp_vec2 _min;
 zp_vec2 _max;
} zp_aabb2d;

/*
 bodies
*/

typedef struct {
 uint16_t         _flags;
 zp_aabb2d        _fit_aabb;
 zp_vec2          _position;
 zp_vec2          _velocity;
 zp_vec2          _force;
 zp_complex       _rotation;
 float            _linear_damping;
 float            _omega;
 float            _torque;
 float            _angular_damping;
 float            _inv_inertia;
 float            _inv_mass;
 float            _restitution;
 float            _friction;
 float            _idle_time;
  
 zp_pool_id       _aabb_node;
 zp_container_id  _id; /* minmize the byte size */
} zp_head2d;

typedef struct {
 zp_head2d _head;
 zp_vec2   _half_size;
} zp_box2d;

typedef union {
 zp_head2d _head;
 zp_box2d  _box;
} zp_body2d; 

zp_cold zp_noinline void zp_body2d_init(zp_body2d *const zp_restrict body, const void *const zp_restrict data);
zp_hot void zp_body2d_updatev(zp_body2d *const zp_restrict body, const void *const zp_restrict world, const float dt);
zp_hot void zp_body2d_updatep(zp_body2d *const zp_restrict body, const void *const zp_restrict world, const float dt);

#endif

