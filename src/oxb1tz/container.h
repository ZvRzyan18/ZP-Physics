#ifndef OXB_CONTAINER_H
#define OXB_CONTAINER_H

#include "oxb1tz/types.h"
#include "oxb1tz/hint.h"

typedef struct {
 oxb_id _val;
} oxb_container_id;

typedef struct {
 oxb_allocator    *_allocator;
 uint8_t          *_bytes;
 oxb_container_id *_id_to_index;
 oxb_container_id *_index_to_id;
 oxb_container_id *_free_list;
 size_t           _size;
 size_t           _max_size;

 size_t           _free_count;
 size_t           _forward_index;

 uint16_t         _stride;
} oxb_container;

extern const oxb_container_id OXB_CONTAINER_NULL_ID;
oxb_inline bool oxb_container_id_isequal(oxb_container_id a, oxb_container_id b) {
 return a._val == b._val;
}
oxb_inline bool oxb_container_id_notequal(oxb_container_id a, oxb_container_id b) {
 return a._val != b._val;
}
oxb_inline bool oxb_container_id_isnull(oxb_container_id a) {
 return oxb_container_id_isequal(a, OXB_CONTAINER_NULL_ID);
}


void oxb_container_create(oxb_container *oxb_restrict const c, oxb_allocator *oxb_restrict const allocator, const uint16_t stride, const size_t reserve);
void oxb_container_destroy(oxb_container *oxb_restrict const c);
oxb_container_id oxb_container_acquire(oxb_container *oxb_restrict const c);
void oxb_container_release(oxb_container *oxb_restrict const c, const oxb_container_id id);
void oxb_container_reset(oxb_container *oxb_restrict const c);

void oxb_container_insertionsort(oxb_container *oxb_restrict const c, bool (*should_swap)(void *a, void * b));


oxb_inline void* oxb_container_get(oxb_container *oxb_restrict const c, const oxb_container_id id) {
 /* catch overflow id */
 assert(id._val < c->_max_size);

 /* catch null */
 assert(!oxb_container_id_isnull(id));
 assert(!oxb_container_id_isnull(c->_id_to_index[id._val]));
 assert(!oxb_container_id_isnull(c->_index_to_id[c->_id_to_index[id._val]._val]));

 return (void*)(c->_bytes + ((size_t)c->_id_to_index[id._val]._val * (size_t)c->_stride));
}

#endif

