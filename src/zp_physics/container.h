#ifndef ZP_CONTAINER_H
#define ZP_CONTAINER_H

#include "zp_physics/types.h"
#include <assert.h>

typedef struct {
 __zp_uid _val;
} zp_container_id;

typedef struct {
 uint8_t  *_bytes;
 zp_container_id *_to_index_lut;
 zp_container_id *_to_id_lut;
 zp_container_id *_free_list;
 size_t          _reserve;
 size_t          _max_size;
 size_t          _size;
 float           _growth_base;
 uint16_t        _stride;
 uint16_t        _increase_count;
} zp_container;

extern const zp_container_id ZP_CONTAINER_NULL_ID;

zp_hot zp_inline int zp_container_id_isequal(zp_container_id a, zp_container_id b) {
 return a._val == b._val;
}
zp_hot zp_inline int zp_container_id_notequal(zp_container_id a, zp_container_id b) {
 return a._val != b._val;
}
zp_hot zp_inline int zp_container_id_isnull(zp_container_id a) {
 return zp_container_id_isequal(a, ZP_CONTAINER_NULL_ID);
}


zp_cold int zp_container_init(zp_container *const zp_restrict c, const uint16_t stride, const size_t reserve, const float growth_base);
zp_cold void zp_container_destroy(zp_container *const zp_restrict c);
zp_container_id zp_container_acquire(zp_container *const zp_restrict c);
void zp_container_release(zp_container *const zp_restrict c, zp_container_id id);
zp_cold void zp_container_insertion_sort(zp_container *const zp_restrict c, int (*should_swap)(void*, void*));

zp_hot zp_inline void* zp_container_get(zp_container *const zp_restrict c, const zp_container_id id) {
 assert(!zp_container_id_isnull(id));
 assert(c->_to_index_lut[id._val]._val < c->_size);
 return (void*)(c->_bytes + ((size_t)c->_to_index_lut[id._val]._val * (size_t)c->_stride));
}


#endif

