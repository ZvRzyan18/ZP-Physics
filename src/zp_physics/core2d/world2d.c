#include "zp_physics/core2d/world2d.h"
#include "zp_physics/core2d/body2d.h"
#include "zp_physics/core2d/solver2d.h"
#include "zp_physics/core2d/collision2d.h"
#include "zp_physics/core2d.h"
#include <string.h>
#include <stdlib.h>
#include <assert.h>

/* condition for body sorting */
static int sort_body_movement(zp_body2d *a, zp_body2d *b) {
	return (a->_head._flags & ZP_BODY_MOVEMENT_MASK_2D) > (b->_head._flags & ZP_BODY_MOVEMENT_MASK_2D);
}

zp_cold int zp_world2d_init(zp_world2d **const zp_restrict world, const zp_create_world2d *const zp_restrict data) {
	*world = (zp_world2d *)malloc(sizeof(zp_world2d));
	zp_world2d *world_instance = *world;

	memset(world_instance, 0, sizeof(zp_world2d));

	world_instance->_dynamic_count = 0;
	world_instance->_kinematic_count = 0;
	world_instance->_sort_body = 0;

	world_instance->_growth_base = data->_growth_base_rate;
	if(zp_unlikely(zp_container_init(&world_instance->_body_container, sizeof(zp_body2d), data->_bodies_initial_reserve, world_instance->_growth_base)))
		return -1;
	if(zp_unlikely(zp_contacthash2d_init(&world_instance->_contacts, data->_hash_max_buckets, data->_contacts_initial_reserve, world_instance->_growth_base)))
		return -1;
	world_instance->_gravity = data->_gravity;
	world_instance->_solver_substeps = data->_solver_substeps;
	world_instance->_time_substeps = data->_time_substeps;

	world_instance->_inv_timestep_substep = 1.0f / (float)world_instance->_time_substeps;
	assert(world_instance->_time_substeps != 0);
	return 0;
}

zp_cold void zp_world2d_destroy(zp_world2d *const zp_restrict world) {
	zp_container_destroy(&world->_body_container);
	zp_contacthash2d_destroy(&world->_contacts);
	free(world);
}

zp_hot void zp_world2d_update(zp_world2d *const zp_restrict world, const float dt) {
	zp_body2d *const start_body = (zp_body2d *)world->_body_container._bytes;

	if(zp_unlikely(world->_sort_body))	{
		zp_container_insertion_sort(&world->_body_container, (void *)sort_body_movement);
		world->_sort_body = 0;
	}

	/* prototype broadphase : bruite force */
	for(uint16_t i = 0; i < world->_body_container._size; i++)	{
		for(uint16_t j = i + 1; j < world->_body_container._size; j++)	{
			zp_body2d *const body_a = start_body + i;
			zp_body2d *const body_b = start_body + j;
			/* early terminate for sleeping object. */
			if((body_a->_head._flags & ZP_BODY_IS_SLEEP_2D) && (body_b->_head._flags & ZP_BODY_IS_SLEEP_2D))
				continue;

			zp_manifold2d m;
			uint8_t intersect = 0;
			memset(&m, 0, sizeof(zp_manifold2d));
			uint8_t queried = zp_contacthash2d_is_queried(&world->_contacts, body_a->_head._id, body_b->_head._id);

			if(queried == 0 || queried == 0xFF) {
				intersect = zp_collision2d_collide(&m, body_a, body_b);
				if(intersect)
					zp_contacthash2d_insert(&world->_contacts, &m);
			}
		}
	}

	zp_contacthash2d_remove_unused(&world->_contacts);

	const float devided_dt = dt * world->_inv_timestep_substep;

	float inv_dt = 1.0f / dt;

	zp_solver_input2d input;
	input._dt = dt;
	input._inv_dt = inv_dt;
	input._solver_substeps = world->_solver_substeps;

	float hertz = 9.33f;
	float damping_ratio = 2.6f;
	float omega = hertz * 6.28318530f;

	/*
  https://box2d.org/posts/2024/02/solver2d/
  effective mass independent soft constraints
 */
	const float a1 = 2.0f * damping_ratio + omega * input._dt;
	const float a2 = input._dt * omega * a1;
	const float a3 = 1.0f / (1.0f + a2);

	input._bias_ratio = omega / a1;
	input._mass_coeff = a2 * a3;
	input._impulse_coeff = a3;

	for(size_t i = 0; i < world->_contacts._memory_pool._size; i++)	{
		zp_contacthash2d_node *const node = ((zp_contacthash2d_node *)world->_contacts._memory_pool._bytes) + i;
		zp_manifold2d *const manifold = &node->_value;
		zp_manifold2d_soft_prepare_contact(manifold, (void *)world, &input);
	}

	/*
  Temporal coherence, the contact data has to be reused
  across several frames.
 */
	for(uint8_t aa = 0; aa < world->_time_substeps; aa++)	{
		for(size_t i = 0; i < world->_contacts._memory_pool._size; i++)	{
			zp_contacthash2d_node *const node = ((zp_contacthash2d_node *)world->_contacts._memory_pool._bytes) + i;
			zp_manifold2d *const manifold = &node->_value;
			zp_manifold2d_soft_presolve_contact(manifold, (void *)world, &input);
		}

		//for(uint16_t i = 0; i < world->_body_container._size; i++) {
		for(size_t i = 0; i < world->_dynamic_count; i++)	{
			zp_body2d *const body_a = start_body + i;
			zp_body2d_updatev(body_a, world, devided_dt);
		}

		for(int substeps = 0; substeps < world->_solver_substeps; substeps++) {
			for(size_t i = 0; i < world->_contacts._memory_pool._size; i++)	{
				zp_contacthash2d_node *const node = ((zp_contacthash2d_node *)world->_contacts._memory_pool._bytes) + i;
				zp_manifold2d *const manifold = &node->_value;
				zp_manifold2d_soft_solve_contact(manifold, (void *)world, &input);
			}
		}

		//for(uint16_t i = 0; i < world->_body_container._size; i++) {
		size_t position_integration_size = world->_dynamic_count + world->_kinematic_count;
		for(size_t i = 0; i < position_integration_size; i++)	{
			zp_body2d *const body_a = start_body + i;
			zp_body2d_updatep(body_a, world, devided_dt);
		}
		
		/* relaxation, improve stability */
		for(size_t i = 0; i < world->_contacts._memory_pool._size; i++)	{
			zp_contacthash2d_node *const node = ((zp_contacthash2d_node *)world->_contacts._memory_pool._bytes) + i;
			zp_manifold2d *const manifold = &node->_value;
			zp_manifold2d_soft_relaxation(manifold, (void *)world, &input);
		}
	}
}

zp_handle zp_world2d_create_body(zp_world2d *const zp_restrict world, const void *const zp_restrict value){
	zp_container_id id = zp_container_acquire(&world->_body_container);
	zp_body2d *const body = (zp_body2d *)zp_container_get(&world->_body_container, id);
	memset(body, 0, sizeof(zp_body2d));
	zp_body2d_init(body, value);
	body->_head._id = id;

	{
		/* sort the body*/
		switch(body->_head._flags & ZP_BODY_MOVEMENT_MASK_2D)	{
		case ZP_BODY_MOVEMENT_DYNAMIC_2D:
			world->_dynamic_count++;
			break;
		case ZP_BODY_MOVEMENT_KINEMATIC_2D:
			world->_kinematic_count++;
			break;
		}
		world->_sort_body = 1;
	}

	return ZP_HANDLE_FLAG_IS_BODY | ((zp_handle)id);
}

void zp_world2d_remove_body(zp_world2d *const zp_restrict world, const zp_handle id)
{
	zp_body2d *const body = (zp_body2d *)zp_container_get(&world->_body_container, (zp_container_id)id);

	{
		switch(body->_head._flags & ZP_BODY_MOVEMENT_MASK_2D) {
		case ZP_BODY_MOVEMENT_DYNAMIC_2D:
			world->_dynamic_count--;
			break;
		case ZP_BODY_MOVEMENT_KINEMATIC_2D:
			world->_kinematic_count--;
			break;
		}
		world->_sort_body = 1;
	}

	zp_container_release(&world->_body_container, (zp_container_id)id);
}

zp_hot void zp_world2d_get_bodydata(zp_world2d *const zp_restrict world, const zp_handle id, zp_bodydata2d *const zp_restrict dat) {
	assert((id & ZP_HANDLE_FLAG_MASK) == ZP_HANDLE_FLAG_IS_BODY);
	zp_body2d *const body = (zp_body2d *)zp_container_get(&world->_body_container, id);
	dat->_position = body->_head._position;

	dat->_rotation = body->_head._rotation;

	dat->_velocity = body->_head._velocity;
	dat->_omega = body->_head._omega;
	dat->_is_sleeping = (body->_head._flags & ZP_BODY_IS_SLEEP_2D) != 0;
}
/*
void insert(ObjectID id, const AABB& box)
{
    int minX = ceil(box.min.x / cellSize);
    int maxX = floor(box.max.x / cellSize);

    int minY = ceil(box.min.y / cellSize);
    int maxY = floor(box.max.y / cellSize);

    for (int y = minY; y <= maxY; ++y)
    {
        for (int x = minX; x <= maxX; ++x)
        {
            grid[hash({x, y})].push_back(id);
        }
    }
}

*/
