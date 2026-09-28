#ifndef ZP_ISLAND2D_H
#define ZP_ISLAND2D_H

#include "zp_physics/core2d/body2d.h"
#include "zp_physics/container.h"
#include "zp_physics/bump.h"

typedef struct {
 size_t          _allocation;
 size_t          _parent;  
 zp_container_id _leader;
 zp_container_id _assistant;
 size_t          _member_count;
 uint16_t        _sleeping_count;
 uint8_t         _invalid;
} zp_island2d_id;

typedef struct {
 zp_container_id _leader;
 size_t          _member_count;
 size_t          _sleeping_count;
} zp_island_representative;

typedef struct {
 zp_container _island_ids;
 zp_bump      _stack;
} zp_island2d;

void zp_island2d_init(zp_island2d *zp_restrict const island, size_t reserve, float growth_base);
void zp_island2d_destroy(zp_island2d *zp_restrict const island);
void zp_island2d_merge(zp_island2d *zp_restrict const island, zp_container *zp_restrict const container, zp_body2d *zp_restrict const body_a, zp_body2d *zp_restrict const body_b);
void zp_island2d_reset(zp_island2d *zp_restrict const island);


#endif

