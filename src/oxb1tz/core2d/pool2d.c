#include "oxb1tz/core2d/pool2d.h"
#include "oxb1tz/core2d/bvh2d.h"
#include <assert.h>

/*
 multi type memory dense pool.
 see @container.c
*/

/*
 OXB_OBJECT_DYNAMIC_2D =   0,
 OXB_OBJECT_KINEMATIC_2D = 1,
 OXB_OBJECT_STATIC_2D =    2,
 OXB_OBJECT_SENSOR_2D =    3,
 OXB_OBJECT_JOINT_2D =     4,
*/
void oxb_pool2d_create(oxb_pool2d *oxb_restrict const pool, oxb_allocator *oxb_restrict const alloc, const size_t initial_reserve) {
 oxb_container_create(&pool->_containers[0], alloc, sizeof(oxb_dynamic2d), initial_reserve);
 oxb_container_create(&pool->_containers[1], alloc, sizeof(oxb_kinematic2d), initial_reserve);
 oxb_container_create(&pool->_containers[2], alloc, sizeof(oxb_static2d), initial_reserve);
}

void oxb_pool2d_destroy(oxb_pool2d *oxb_restrict const pool) {
 oxb_container_destroy(&pool->_containers[0]);
 oxb_container_destroy(&pool->_containers[1]);
 oxb_container_destroy(&pool->_containers[2]);
}


oxb_func oxb_handle oxb_pool2d_acquire(oxb_pool2d *oxb_restrict const pool, const oxb_object2d_info *oxb_restrict const info) {
 size_t pool_index = (info->_flags & OXB_OBJECT_MASK_2D) >> 26;
 assert(pool_index < 3);
 oxb_container_id container_id = oxb_container_acquire(&pool->_containers[pool_index]);
 oxb_handle handle;
 handle._type = pool_index;
 handle._id = container_id._val;
 
 oxb_bvh2d *bvh = (oxb_bvh2d*)pool->_bvh;

 switch(handle._type) {
  case 0:
  {
  	oxb_dynamic2d *dynamic_body = (oxb_dynamic2d*)oxb_container_get(&pool->_containers[0], container_id);
  	oxb_dynamic2d_init(dynamic_body, info);
  	dynamic_body->_allocation = handle;
  	
  	dynamic_body->_bvh_id = oxb_bvh2d_insert(bvh, dynamic_body->_aabb, handle);
  }
  break;
  case 1:
  {
  	oxb_kinematic2d *kinematic_body = (oxb_kinematic2d*)oxb_container_get(&pool->_containers[1], container_id);
  	oxb_kinematic2d_init(kinematic_body, info);
  	kinematic_body->_allocation = handle;
  	
  	kinematic_body->_bvh_id = oxb_bvh2d_insert(bvh, kinematic_body->_aabb, handle);
  }
  break;
  case 2:
  {
  	oxb_static2d *static_body = (oxb_static2d*)oxb_container_get(&pool->_containers[2], container_id);
  	oxb_static2d_init(static_body, info);
  	static_body->_allocation = handle;
  	
  	static_body->_bvh_id = oxb_bvh2d_insert(bvh, static_body->_aabb, handle);
  }
  break;
  default:
   assert(0);
   oxb_unreachable();
 }
 return handle;
}


oxb_func void oxb_pool2d_release(oxb_pool2d *oxb_restrict const pool, const oxb_handle id) {
 assert(id._id < 5);
 oxb_container_id container_id;
 container_id._val = id._id;
 oxb_container_release(&pool->_containers[id._type], container_id);

 oxb_bvh2d *bvh = (oxb_bvh2d*)pool->_bvh;

 switch(id._type) {
  case 0:
  {
  	oxb_dynamic2d *dynamic_body = (oxb_dynamic2d*)oxb_container_get(&pool->_containers[0], container_id);
  	oxb_bvh2d_remove(bvh, dynamic_body->_bvh_id);
  }
  break;
  case 1:
  {
  	oxb_kinematic2d *kinematic_body = (oxb_kinematic2d*)oxb_container_get(&pool->_containers[1], container_id);
  	oxb_bvh2d_remove(bvh, kinematic_body->_bvh_id);
  }
  break;
  case 2:
  {
  	oxb_static2d *static_body = (oxb_static2d*)oxb_container_get(&pool->_containers[2], container_id);
  	oxb_bvh2d_remove(bvh, static_body->_bvh_id);
  }
  break;
 }
}


oxb_func void oxb_pool2d_getdata(oxb_pool2d *oxb_restrict const pool, oxb_object2d_data *oxb_restrict const dat, const oxb_handle id) {
 switch(id._type) {
  case 0:
  {
  	oxb_container_id lookup;
  	lookup._val = id._id;
  	oxb_dynamic2d *dynamic_body = oxb_container_get(&pool->_containers[0], lookup);
  	
  	dat->_position = oxb_vec2_load(dynamic_body->_position._arr);
  	dat->_rotation = oxb_vec2_load(dynamic_body->_rotation._arr);
  	dat->_velocity = oxb_vec2_load(dynamic_body->_velocity._arr);
  	dat->_omega = dynamic_body->_omega;
  }
  break;
  case 1:
  break;
  case 2:
   assert(0);
   oxb_unreachable();
  break;
  default:
   assert(0);
   oxb_unreachable();
 }
}

