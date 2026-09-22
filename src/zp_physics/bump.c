#include "zp_physics/bump.h"
#include "zp_physics/math.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>

zp_hot zp_noinline static int bump_resize(zp_bump *const zp_restrict b) {
 zp_bump new_bump;

 new_bump._max_size = b->_max_size + (uint16_t)zp_ceil(zp_exp2(((float)++b->_increase_count) * b->_growth_base));
 new_bump._growth_base = b->_growth_base;
 new_bump._increase_count = b->_increase_count;
 
 new_bump._bytes = (uint8_t*)malloc(new_bump._max_size);
 if(zp_unlikely(!new_bump._bytes))
  return -1;
  
 new_bump._size = b->_size;
 memcpy(new_bump._bytes, b->_bytes, b->_size);
 zp_bump_destroy(b);
 memcpy(b, &new_bump, sizeof(zp_bump));
 return 0;
}

/*
 bump/stack-based allocator
*/
zp_cold int zp_bump_init(zp_bump *const zp_restrict b, const size_t initial_size, const float growth_base) {
 b->_bytes = (uint8_t*)malloc(initial_size);
 b->_max_size = initial_size;
 b->_size = 0;
 
 b->_growth_base = zp_log2(growth_base);
 b->_increase_count = (uint16_t)zp_ceil(zp_log2((float)initial_size) / b->_growth_base);

 if(zp_unlikely(!b->_bytes))
  return -1;
 return 0;
}


zp_cold void zp_bump_destroy(zp_bump *const zp_restrict b) {
 free(b->_bytes);
}


size_t zp_bump_acquire(zp_bump *const zp_restrict b, const uint16_t requested_bytes, const uint8_t alignment) {
 size_t current_pointer;
 current_pointer = zp_alignto(b->_size, alignment);
 b->_size = current_pointer + requested_bytes;
 if(b->_size >= b->_max_size) {
  if(bump_resize(b))
   return 0xFFFFFFFF;
 }
 return current_pointer;
}



void zp_bump_reset(zp_bump *const zp_restrict b, const size_t offset) {
 b->_size = offset;
 assert(b->_size < b->_max_size);
}



