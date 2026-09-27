#ifndef ZP_ISLAND2D_H
#define ZP_ISLAND2D_H

#include "zp_physics/core2d/body2d.h"
#include "zp_physics/container.h"

typedef struct {
 zp_container_id _leader;
 uint16_t        _member_count;
 uint16_t        _sleeping_count;
 uint8_t         _invalid;
} zp_island2d_id;

typedef struct {
 zp_container _island_ids;
} zp_island2d;

void zp_island2d_init(zp_island2d *zp_restrict const island, size_t reserve, float growth_base);
void zp_island2d_destroy(zp_island2d *zp_restrict const island);


#endif

