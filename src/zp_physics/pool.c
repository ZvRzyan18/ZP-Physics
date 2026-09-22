#include "zp_physics/pool.h"
#include "zp_physics/math.h"
#include <string.h>
#include <stdlib.h>

static_assert(sizeof(zp_pool_id) == ZP_ID_SIZE, "not expected size");
const zp_pool_id ZP_POOL_NULL_ID = {._val = ZP_ID_MAX};


zp_noinline zp_cold static int resize(zp_pool *const zp_restrict pool) {
 zp_pool new_pool;
 
 /* not exact, but still it can grow in a decent way, and able to minimize the realloc */
 new_pool._reserve = pool->_reserve + (uint16_t)zp_ceil(zp_exp2(((float)++pool->_increase_count) * pool->_growth_base));
 assert(new_pool._reserve > pool->_reserve);

 new_pool._max_size = pool->_max_size + new_pool._reserve;
 new_pool._size = pool->_size;
 new_pool._stride = pool->_stride;
 new_pool._increase_count = pool->_increase_count;
 new_pool._growth_base = pool->_growth_base;
 
 const size_t free_list_size = (size_t)new_pool._max_size * sizeof(zp_pool_id);
 const size_t bytes_size = (size_t)new_pool._max_size * (size_t)new_pool._stride;
 new_pool._bytes = (uint8_t*)malloc(free_list_size + bytes_size + ZP_MEMORY_ALIGNMENT);
 if(zp_unlikely(!new_pool._bytes))
  return -1;
 new_pool._free_list = (zp_pool_id*)zp_alignto((uint64_t)(new_pool._bytes + bytes_size), ZP_MEMORY_ALIGNMENT);

 assert((((uint64_t)new_pool._bytes) % 2) == 0);
 assert((((uint64_t)new_pool._free_list) % ZP_MEMORY_ALIGNMENT) == 0);

 memcpy(new_pool._bytes, pool->_bytes, pool->_max_size * pool->_stride);
 memcpy(new_pool._free_list, pool->_free_list, pool->_max_size * sizeof(zp_pool_id));

 for(size_t i = pool->_max_size; i < new_pool._max_size; i++)
  new_pool._free_list[i]._val = i;
 
 zp_pool_destroy(pool);
 memcpy(pool, &new_pool, sizeof(zp_pool));
 return 0;
}


zp_cold int zp_pool_init(zp_pool *const zp_restrict pool, const uint16_t stride, const uint16_t reserve, const float growth_base) {
 assert(reserve != 0);
 assert(stride != 0);

 /* not exact, but still it can grow in a decent way, and able to minimize the realloc */
 pool->_growth_base = zp_log2(growth_base);
 pool->_increase_count = (uint16_t)zp_ceil(zp_log2((float)reserve) / pool->_growth_base);
 
 pool->_reserve = reserve;
 pool->_max_size = reserve;
 pool->_size = 0;
 pool->_stride = stride;
 
 const size_t free_list_size = pool->_max_size * sizeof(zp_pool_id);
 const size_t bytes_size = pool->_max_size * stride;
 pool->_bytes = (uint8_t*)malloc(free_list_size + bytes_size + ZP_MEMORY_ALIGNMENT);
 if(zp_unlikely(!pool->_bytes))
  return -1;
 pool->_free_list = (zp_pool_id*)zp_alignto((uint64_t)(pool->_bytes) + bytes_size, ZP_MEMORY_ALIGNMENT);

 assert((((uint64_t)pool->_bytes) % 2) == 0);
 assert((((uint64_t)pool->_free_list) % ZP_MEMORY_ALIGNMENT) == 0);

 for(size_t i = 0; i < pool->_max_size; i++)
  pool->_free_list[i]._val = i;
  
 return 0;
}


zp_cold void zp_pool_destroy(zp_pool *const zp_restrict pool) {
 free(pool->_bytes);
}



zp_pool_id zp_pool_acquire(zp_pool *const zp_restrict pool) {
 if(pool->_size >= pool->_max_size)
  if(zp_unlikely(resize(pool)))
   return ZP_POOL_NULL_ID;
 return pool->_free_list[pool->_size++];
}


void zp_pool_release(zp_pool *const zp_restrict pool, const zp_pool_id id) {
 assert(pool->_size != 0);
 assert(!zp_pool_id_isnull(id));
 pool->_free_list[--pool->_size] = id;
}



