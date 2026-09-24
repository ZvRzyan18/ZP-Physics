#ifndef ZP_BROADPHASE_H
#define ZP_BROADPHASE_H

#include "zp_physics/vector.h"
#include "zp_physics/pool.h"
#include "zp_physics/container.h"
#include "zp_physics/container.h"
#include "zp_physics/bump.h"
#include "zp_physics/core2d/body2d.h"


typedef struct {
 zp_aabb2d _aabb;
 zp_pool_id _left, _right, _parent, _id;
 zp_container_id _body;
} zp_broadphase2d_node;

typedef struct {
 uint32_t _node_visits;
 uint32_t _max_depth;
 uint32_t _pair_test;
} zp_broadphase2d_treetraversalinfo;

typedef struct {
 zp_pool    _node_allocator;
 zp_bump    _stack;
 float      _aabb_margin;
 uint32_t   _leaf_size;
 zp_pool_id _root;
 uint8_t    _frame_tolerance;
} zp_broadphase2d;

int zp_broadphase2d_init(zp_broadphase2d *zp_restrict const bp, const float growth_base, const float aabb_margin, const size_t stack_initial_size);
void zp_broadphase2d_destroy(zp_broadphase2d *zp_restrict const bp);
void zp_broadphase2d_update_element(zp_broadphase2d *zp_restrict const bp, const zp_pool_id id, const zp_aabb2d fit);
zp_pool_id zp_broadphase2d_insert_element(zp_broadphase2d *zp_restrict const bp, const zp_aabb2d fit, const zp_container_id body_id);
void zp_broadphase2d_remove_element(zp_broadphase2d *zp_restrict const bp, const zp_pool_id id);

void zp_broadphase2d_traverse_pairs(zp_broadphase2d *zp_restrict const bp, void (*func)(zp_container_id, zp_container_id, void*), void *ptr);
void zp_broadphase2d_greedy_rebuild_tree(zp_broadphase2d *zp_restrict const bp, void *zp_restrict w);
void zp_broadphase2d_optimize_tree(zp_broadphase2d *zp_restrict const bp, void *zp_restrict w);
int zp_broadphase2d_should_rebuild(zp_broadphase2d *zp_restrict const bp, uint8_t frame_tolerance);

#endif


