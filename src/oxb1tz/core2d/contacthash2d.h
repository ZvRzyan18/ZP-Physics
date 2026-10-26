#ifndef OXB_CONTACTHASH2D_H
#define OXB_CONTACTHASH2D_H

#include "oxb1tz/core2d/manifold2d.h"
#include "oxb1tz/container.h"

typedef struct {
 oxb_manifold2d _manifold;
 oxb_container_id _next;
 oxb_container_id _prev;
 oxb_container_id _allocation;
} oxb_contacthash2d_node;


typedef struct {
 oxb_container    _pool;
 oxb_container_id *_buckets;
 size_t           _buckets_size;
 size_t           _buckets_mask;
} oxb_contacthash2d;


void oxb_contacthash2d_create(oxb_contacthash2d *oxb_restrict const c, oxb_allocator *oxb_restrict const alloc, const size_t reserve);
void oxb_contacthash2d_destroy(oxb_contacthash2d *oxb_restrict const c);
void oxb_contacthash2d_insert(oxb_contacthash2d *oxb_restrict const c, const oxb_manifold2d *oxb_restrict const m);
oxb_contacthash2d_node* oxb_contacthash2d_get(oxb_contacthash2d *oxb_restrict const c, const oxb_handle a, const oxb_handle b);
void oxb_contacthash2d_remove(oxb_contacthash2d *oxb_restrict const c, const oxb_contacthash2d_node *oxb_restrict const n);

#endif


