#ifndef OXB_ISLAND2D_H
#define OXB_ISLAND2D_H

#include "oxb1tz/types.h"
#include "oxb1tz/core2d/pool2d.h"

oxb_handle oxb_island2d_find(oxb_pool2d *oxb_restrict const pool, const oxb_handle h);
void oxb_island2d_merge(oxb_pool2d *oxb_restrict const pool, const oxb_handle a, const oxb_handle b);

#endif

