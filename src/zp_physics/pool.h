#ifndef ZP_POOL_H
#define ZP_POOL_H

#include "zp_physics/types.h"
#include <assert.h>


typedef struct {
 __zp_uid _val;
} zp_pool_id;

typedef struct {
 uint8_t    *_bytes;
 zp_pool_id *_free_list;
 size_t      _size;
 size_t      _reserve;
 size_t      _max_size;
 float       _growth_base;
 uint16_t    _stride;
 uint16_t    _increase_count;
} zp_pool;


extern const zp_pool_id ZP_POOL_NULL_ID;

zp_hot zp_inline int zp_pool_id_isequal(zp_pool_id a, zp_pool_id b) {
 return a._val == b._val;
}
zp_hot zp_inline int zp_pool_id_notequal(zp_pool_id a, zp_pool_id b) {
 return a._val != b._val;
}
zp_hot zp_inline int zp_pool_id_isnull(zp_pool_id id) {
 return zp_pool_id_isequal(id, ZP_POOL_NULL_ID);
}


zp_cold int zp_pool_init(zp_pool *const zp_restrict pool, const uint16_t stride, const uint16_t reserve, const float growth_base);
zp_cold void zp_pool_destroy(zp_pool *const zp_restrict pool);
zp_pool_id zp_pool_acquire(zp_pool *const zp_restrict pool);
void zp_pool_release(zp_pool *const zp_restrict pool, const zp_pool_id id);
zp_inline void* zp_pool_get(zp_pool *const zp_restrict pool, zp_pool_id id) {
 assert(!zp_pool_id_isnull(id));
 return pool->_bytes + ((size_t)id._val) * ((size_t)pool->_stride);
}


#endif

