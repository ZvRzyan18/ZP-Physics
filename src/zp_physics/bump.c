#include "zp_physics/bump.h"
#include "zp_physics/math.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>

zp_hot zp_noinline static void bump_resize(zp_bump *zp_restrict b) {
 zp_bump new_bump;
 new_bump._max_size = b->_size * 2;
 new_bump._bytes = (uint8_t*)malloc(new_bump._max_size);
 new_bump._size = b->_size;
 memcpy(new_bump._bytes, b->_bytes, b->_size);
 zp_bump_destroy(b);
 memcpy(b, &new_bump, sizeof(zp_bump));
}

/*
 bump/stack-based allocator
*/
zp_cold int zp_bump_init(zp_bump *zp_restrict b, size_t initial_size) {
 b->_bytes = (uint8_t*)malloc(initial_size);
 b->_max_size = initial_size;
 b->_size = 0;
 
 if(zp_unlikely(!b->_bytes))
  return -1;
 return 0;
}


zp_cold void zp_bump_destroy(zp_bump *zp_restrict b) {
 free(b->_bytes);
}


size_t zp_bump_acquire(zp_bump *zp_restrict b, uint16_t requested_bytes, uint8_t alignment) {
 size_t current_pointer;
 current_pointer = zp_alignto(b->_size, alignment);
 b->_size = current_pointer + requested_bytes;
 if(b->_size >= b->_max_size) {
  bump_resize(b);
 }
 return current_pointer;
}



void zp_bump_reset(zp_bump *zp_restrict b, size_t offset) {
 b->_size = offset;
 assert(b->_size < b->_max_size);
}



