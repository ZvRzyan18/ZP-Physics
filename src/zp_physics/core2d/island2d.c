#include "zp_physics/core2d/island2d.h"
#include <assert.h>
#include <string.h>

/**********************************************
*               ISLAND OPERATION              *
***********************************************/




void zp_island2d_init(zp_island2d *zp_restrict const island, size_t reserve, float growth_base) {
 zp_container_init(&island->_island_ids, sizeof(zp_island2d_id), reserve, growth_base);
}

void zp_island2d_destroy(zp_island2d *zp_restrict const island) {
 zp_container_destroy(&island->_island_ids);
}



