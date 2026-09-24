#ifndef ZP_CORE2D_H
#define ZP_CORE2D_H

#include "zp_physics/types.h"
#include "zp_physics/hint.h"

/* body states */
#define ZP_BODY_MOVEMENT_MASK_2D      0b1110000000000000
#define ZP_BODY_MOVEMENT_DYNAMIC_2D   0b0010000000000000
#define ZP_BODY_MOVEMENT_KINEMATIC_2D 0b0100000000000000
#define ZP_BODY_MOVEMENT_STATIC_2D    0b0110000000000000

#define ZP_BODY_IS_SLEEP_2D           0b0001000000000000
#define ZP_BODY_LOCK_POS_X_2D         0b0000100000000000
#define ZP_BODY_LOCK_POS_Y_2D         0b0000010000000000
#define ZP_BODY_LOCK_ROT_2D           0b0000001000000000

#define ZP_BODY_MASK_2D               0b0000000000111111
#define ZP_BODY_BOX_2D                0b0000000000000001
#define ZP_BODY_CIRCLE_2D             0b0000000000000010

#define ZP_JOINT_MASK_2D              0b0000000000111111
#define ZP_JOINT_DISTANCE_2D          0b0000000000000001
#define ZP_JOINT_FIXED_2D             0b0000000000000010

/* external handle states. */
#define ZP_HANDLE_FLAG_IS_BODY_2D     0b00000001
#define ZP_HANDLE_FLAG_IS_JOINT_2D    0b00000010
#define ZP_HANDLE_FLAG_IS_SENSOR_2D   0b00000011

#define ZP_WORLD_DT_LIMIT_CLAMP_2D     1
#define ZP_WORLD_DT_LIMIT_ADAPTIVE_2D  2
#define ZP_WORLD_DT_LIMIT_NONE_2D      3


/*
 create infos
*/


typedef struct {
 struct {
  size_t    _broadphase_stack_bytes_inital_reserve;
  float     _growth_base_rate;
  uint16_t  _bodies_initial_reserve;
  uint16_t  _joints_initial_reserve;
  uint16_t  _contacts_initial_reserve;
  uint16_t  _max_buckets_initial_reserve;
 } _memory;
 
 struct {
  float   _aabb_margin;
  uint8_t _frame_tolerance;
 } _broadphase_info;
 
 struct {
  float   _hertz;
  float   _damping_ratio;
 } _contact_constraint;
 
 
 zp_vec2  _gravity;
 uint8_t  _dt_limit;
 uint8_t  _solver_substeps;
 uint8_t  _time_substeps;
} zp_create_world2d;


typedef struct {
 uint16_t _flags;
 zp_vec2  _position;
 zp_vec2  _velocity;
 zp_vec2  _force;
 zp_vec2  _half_size;
 float    _linear_damping;
 float    _rotation;
 float    _omega;
 float    _torque;
 float    _density;
 float    _angular_damping;
 float    _restitution;
 float    _friction;
} zp_create_body2d;



typedef struct {
 zp_vec2    _position;
 zp_complex _rotation;
 zp_vec2    _velocity;
 float      _omega;
 uint8_t    _is_sleeping;
} zp_bodydata2d;

typedef struct zp_world2d zp_world2d;

#ifdef __cplusplus
extern "C" {
#endif

zp_cold int zp_world2d_init(zp_world2d **const zp_restrict world, const zp_create_world2d *const zp_restrict data);
zp_cold void zp_world2d_destroy(zp_world2d *const zp_restrict world);
zp_hot void zp_world2d_update(zp_world2d *const zp_restrict world, const float dt);
zp_handle zp_world2d_create_body(zp_world2d *const zp_restrict world, const void *const zp_restrict value);
void zp_world2d_remove_body(zp_world2d *const zp_restrict world, const zp_handle id);
zp_hot void zp_world2d_get_bodydata(zp_world2d *const zp_restrict world, const zp_handle id, zp_bodydata2d *const zp_restrict dat);
zp_cold void zp_world2d_greedy_rebuild_tree(zp_world2d *const zp_restrict world);
zp_cold void zp_world2d_optimize_rebuild_tree(zp_world2d *const zp_restrict world);
zp_cold void zp_world2d_optimize_body_container(zp_world2d *const zp_restrict world);

#ifdef __cplusplus
}
#endif


#endif


