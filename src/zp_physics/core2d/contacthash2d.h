#ifndef ZP_CONTACTHASH2D_H
#define ZP_CONTACTHASH2D_H

#include "zp_physics/hint.h"
#include "zp_physics/container.h"
#include "zp_physics/core2d/manifold2d.h"

typedef struct {
 zp_manifold2d _value;
 zp_container_id _allocation;
 zp_container_id _next;
 zp_container_id _prev;
 uint8_t  _queried;
} zp_contacthash2d_node;

typedef struct {
 zp_container     _memory_pool;
 zp_container_id  *_bucket;
 size_t           _bucket_index_mask;
 size_t           _bucket_size;
} zp_contacthash2d;



zp_cold int zp_contacthash2d_init(zp_contacthash2d *const zp_restrict hash, const size_t bucket_size, const size_t reserve, const float growth_base);
zp_cold void zp_contacthash2d_destroy(zp_contacthash2d *const zp_restrict hash);
void zp_contacthash2d_insert(zp_contacthash2d *const zp_restrict hash, const zp_manifold2d *const zp_restrict m);
void zp_contacthash2d_remove_unused(zp_contacthash2d *const zp_restrict hash);
uint8_t zp_contacthash2d_is_queried(zp_contacthash2d *const zp_restrict hash, const zp_container_id a, const zp_container_id b);


#endif


