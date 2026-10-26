#ifndef OXB_POOL2D_H
#define OXB_POOL2D_H

#include "oxb1tz/container.h"
#include "oxb1tz/core2d/object2d.h"
#include "oxb1tz/core2d.h"

struct oxb_pool2d {
 oxb_container _containers[5];
 void *_bvh;
};

void oxb_pool2d_create(oxb_pool2d *oxb_restrict const pool, oxb_allocator *oxb_restrict const alloc, const size_t initial_reserve);
void oxb_pool2d_destroy(oxb_pool2d *oxb_restrict const pool);

oxb_func oxb_handle oxb_pool2d_acquire(oxb_pool2d *oxb_restrict const pool, const oxb_object2d_info *oxb_restrict const info);
oxb_func void oxb_pool2d_release(oxb_pool2d *oxb_restrict const pool, const oxb_handle id);
oxb_func void oxb_pool2d_getdata(oxb_pool2d *oxb_restrict const pool, oxb_object2d_data *oxb_restrict const dat, const oxb_handle id);

oxb_inline void* oxb_pool2d_get(oxb_pool2d *oxb_restrict const pool, const oxb_handle id) {
 assert(id._type < 3);
 oxb_container_id container_id;
 container_id._val = id._id;
 return oxb_container_get(&pool->_containers[id._type], container_id);
}

oxb_inline void oxb_pool2d_collider(oxb_pool2d *oxb_restrict const pool, oxb_collider2d *oxb_restrict const coll, const oxb_handle id) {
 oxb_container_id lookup;
 lookup._val = id._id;
 switch(id._type) {
  case 0:
   oxb_dynamic2d_collider(coll, (oxb_dynamic2d*)oxb_container_get(&pool->_containers[0], lookup));
  break;
  case 1:
   oxb_kinematic2d_collider(coll, (oxb_kinematic2d*)oxb_container_get(&pool->_containers[1], lookup));
  break;
  case 2:
   oxb_static2d_collider(coll, (oxb_static2d*)oxb_container_get(&pool->_containers[2], lookup));
  break;
  default:
   assert(0);
   oxb_unreachable();
 }
}

#endif
