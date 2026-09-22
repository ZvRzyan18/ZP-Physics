#ifndef ZP_COLLISION2D_H
#define ZP_COLLISION2D_H

#include "zp_physics/hint.h"
#include "zp_physics/core2d/manifold2d.h"

uint8_t zp_collision2d_collide(zp_manifold2d *const zp_restrict out, const zp_body2d *zp_restrict a, const zp_body2d *zp_restrict b);

#endif

