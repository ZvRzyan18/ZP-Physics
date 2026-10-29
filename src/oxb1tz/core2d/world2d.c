#include "oxb1tz/core2d/world2d.h"
#include "oxb1tz/memory.h"



oxb_noinline void broadphase_query(oxb_handle ha, oxb_handle hb, void *w) {
 oxb_world2d *world = (oxb_world2d*)w;

 oxb_collider2d a, b;
 oxb_pool2d_collider(&world->_pool2d, &a, ha);
 oxb_pool2d_collider(&world->_pool2d, &b, hb);

 if((a._aabb._max._arr[0] < b._aabb._min._arr[0] || a._aabb._min._arr[0] > b._aabb._max._arr[0] || a._aabb._max._arr[1] < b._aabb._min._arr[1] || a._aabb._min._arr[1] > b._aabb._max._arr[1]))
  return;


 oxb_contacthash2d_node *node = oxb_contacthash2d_get(&world->_contacts, ha, hb);
 
 /* skip duplicates */
 if(node != NULL) {
  if(node->_manifold._flags & OXB_MANIFOLD2D_QUERIED) {
   return;
  }
 }
 
 /* narrowphase */
 oxb_manifold2d m;
 oxb_memory_set(&m, 0x00, sizeof(oxb_manifold2d));

 bool result = oxb_manifold2d_object_vs_object[a._flags][b._flags](&m, &a, &b);
 if(!result)
  return;
 
 if(result) {
  if(node == NULL) {
   m._flags |= OXB_MANIFOLD2D_QUERIED;
   oxb_contacthash2d_insert(&world->_contacts, &m);
  } else {
   node->_manifold._flags |= OXB_MANIFOLD2D_QUERIED;
   oxb_manifold2d_combine(&node->_manifold, &m);
  }
 }
}


oxb_world2d* oxb_world2d_create(const oxb_world2d_info *oxb_restrict const info) {
 oxb_world2d *world = (oxb_world2d*)oxb_aligned_alloc((oxb_allocator*)&info->_allocator, sizeof(oxb_world2d));
 if(oxb_unlikely(!world)) {
  OXB_PANIC("oxb_world2d_create() : allocation failed. ");
 }
 oxb_memory_set(world, 0, sizeof(oxb_world2d));

 world->_allocator = info->_allocator;
 world->_gravity = info->_gravity;
 world->_substeps = info->_substeps;
 oxb_pool2d_create(&world->_pool2d, &world->_allocator, 2);
 oxb_contacthash2d_create(&world->_contacts, &world->_allocator, 32);
 oxb_bvh2d_create(&world->_bvh2d, &world->_allocator, 32, 1024, 2.0f);

 world->_pool2d._world = (void*)world;
 return world;
}


void oxb_world2d_destroy(oxb_world2d *oxb_restrict const world) {
 oxb_bvh2d_destroy(&world->_bvh2d);
 oxb_contacthash2d_destroy(&world->_contacts);
 oxb_pool2d_destroy(&world->_pool2d);
 oxb_aligned_free(&world->_allocator, world);
}


oxb_pool2d* oxb_world2d_getpool2d(oxb_world2d *oxb_restrict const world) {
 return &world->_pool2d;
}


void oxb_world2d_step(oxb_world2d *oxb_restrict const world, const float dt) {
 world->_dt = dt;
 float divided_dt = dt / (float)world->_substeps;

 if(oxb_bvh2d_shouldrebuild(&world->_bvh2d, 5)) 
  oxb_bvh2d_mediansplit(&world->_bvh2d, &world->_pool2d);

 oxb_bvh2d_traversepairs(&world->_bvh2d, broadphase_query, world);


 float damping_ratio = 2.6f;
 float hertz = 9.33f;

 float soft_dt = dt;
 
 float omega = hertz * 6.28318530f;
  
	const float a1 = 2.0f * damping_ratio + omega * soft_dt;
	const float a2 = soft_dt * omega * a1;
	const float a3 = 1.0f / (1.0f + a2);
	world->_input._bias_coeff = omega / a1;
	world->_input._mass_coeff = a2 * a3;
	world->_input._impulse_coeff = a3;
	world->_input._inv_dt = 1.0f / soft_dt;
	world->_input._ctx = (void*)world;


 
 /* remove unused contacts. */
 for(size_t i = 0; i < world->_contacts._pool._size; ) {
 	oxb_contacthash2d_node *node = ((oxb_contacthash2d_node*)world->_contacts._pool._bytes) + i;
  if(node->_manifold._flags & OXB_MANIFOLD2D_QUERIED) {
   node->_manifold._flags = 0;
   i++;
  } else {
   oxb_contacthash2d_remove(&world->_contacts, node);
  }
 }

 oxb_memory_set(world->_colors, 0, sizeof(uint32_t) * 32);
 world->_solo_manifolds = 0;
 for(size_t i = 0; i < world->_contacts._pool._size; i++) {
 	oxb_contacthash2d_node *node = ((oxb_contacthash2d_node*)world->_contacts._pool._bytes) + i;
  oxb_solver2d_contact_prepare(&node->_manifold, &world->_pool2d, &world->_input);
 }



 for(size_t iteration = 0; iteration < world->_substeps; iteration++) {

 for(size_t i = 0; i < world->_contacts._pool._size; i++) {
 	oxb_contacthash2d_node *node = ((oxb_contacthash2d_node*)world->_contacts._pool._bytes) + i;
  oxb_solver2d_contact_presolve(&node->_manifold, &world->_pool2d, &world->_input);
 }

 for(size_t i = 0; i < world->_contacts._pool._size; i++) {
 	oxb_contacthash2d_node *node = ((oxb_contacthash2d_node*)world->_contacts._pool._bytes) + i;
  oxb_solver2d_contact_solve(&node->_manifold, &world->_pool2d, &world->_input);
 }


 for(size_t i = 0; i < world->_pool2d._containers[0]._size; ++i) {
  oxb_dynamic2d *dynamic_body = ((oxb_dynamic2d*)world->_pool2d._containers[0]._bytes) + i;
  oxb_dynamic2d_integrate(dynamic_body, world->_gravity, divided_dt);
  oxb_bvh2d_update(&world->_bvh2d, dynamic_body->_bvh_id, dynamic_body->_aabb);
 }

 for(size_t i = 0; i < world->_contacts._pool._size; i++) {
 	oxb_contacthash2d_node *node = ((oxb_contacthash2d_node*)world->_contacts._pool._bytes) + i;
  oxb_solver2d_contact_solve(&node->_manifold, &world->_pool2d, &world->_input);
 }
 
 } //iteration

}

