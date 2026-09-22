#include "zp_physics/container.h"
#include "zp_physics/math.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

static_assert(sizeof(zp_container_id) == ZP_ID_SIZE, "not expected size");
const zp_container_id ZP_CONTAINER_NULL_ID = {._val = ZP_ID_MAX};

zp_hot zp_inline void memswap(void *zp_restrict const a, void *zp_restrict const b, const size_t stride) {
 assert(stride < 512);
 
 uint8_t tmp_buff[512]; /* more portable than vla */
 memcpy(tmp_buff, a, stride);
 memcpy(a, b, stride);
 memcpy(b, tmp_buff, stride);
}


zp_cold zp_noinline static int resize(zp_container *const zp_restrict c) {
 zp_container new_container;
 /* not exact, but still it can grow in a decent way, and able to minimize the realloc */

	new_container._reserve = c->_reserve + (size_t)zp_ceil(zp_exp2(((float)++c->_increase_count) * c->_growth_base));
 assert(new_container._reserve > c->_reserve);
		
	new_container._max_size = c->_size + new_container._reserve;
	new_container._stride = c->_stride;
	new_container._size = c->_size;
	
	new_container._increase_count = c->_increase_count;
 new_container._growth_base = c->_growth_base;

	const size_t bytes_size = new_container._max_size * new_container._stride;
	const size_t helper_size = sizeof(zp_container_id) * new_container._max_size;
	
	const size_t luts_size = helper_size * 3 + bytes_size;

 uint8_t *ptr = (uint8_t*)malloc(luts_size + ZP_MEMORY_ALIGNMENT * 3);
 new_container._bytes = ptr;
 ptr += bytes_size;
 new_container._to_index_lut = (zp_container_id*)zp_alignto((uint64_t)ptr, ZP_MEMORY_ALIGNMENT);
 ptr = ((uint8_t*)new_container._to_index_lut) + helper_size;
 new_container._to_id_lut = (zp_container_id*)zp_alignto((uint64_t)ptr, ZP_MEMORY_ALIGNMENT);
 ptr = ((uint8_t*)new_container._to_id_lut) + helper_size;
 new_container._free_list = (zp_container_id*)zp_alignto((uint64_t)ptr, ZP_MEMORY_ALIGNMENT);

 assert(((uint64_t)ptr - (uint64_t)new_container._bytes) < (luts_size + ZP_MEMORY_ALIGNMENT * 3) );
 assert((((uint64_t)new_container._bytes) % 2) == 0);
 assert((((uint64_t)new_container._to_index_lut) % ZP_MEMORY_ALIGNMENT) == 0);
 assert((((uint64_t)new_container._to_id_lut) % ZP_MEMORY_ALIGNMENT) == 0);
 assert((((uint64_t)new_container._free_list) % ZP_MEMORY_ALIGNMENT) == 0);

 if(zp_unlikely(!new_container._bytes))
  return -1;
 
	const size_t _bytes_size = c->_max_size * new_container._stride;
	const size_t _helper_size = sizeof(zp_container_id) * c->_max_size;

 memcpy(new_container._bytes, c->_bytes, _bytes_size);
 memcpy(new_container._free_list, c->_free_list, _helper_size);
 memcpy(new_container._to_id_lut, c->_to_id_lut, _helper_size);
 memcpy(new_container._to_index_lut, c->_to_index_lut, _helper_size);

 size_t i = c->_max_size;
 while(i < new_container._max_size) {
  new_container._to_index_lut[i]._val = i;
  new_container._free_list[i]._val = i;
  i++;
 }
 zp_container_destroy(c);
 memcpy(c, &new_container, sizeof(zp_container));
 return 0;
}


/*

 lets say we have 
 
 Elements [A, E]
 Indices  [0, 4]
 
 array         ids
 •---•  •---•
 | A |  | 0 |   0
 •---•  •---•
 | B |  | 1 |   1
 •---•  •---•
 | C |  | 2 |   2
 •---•  •---•
 | D |  | 3 |   3
 •---•  •---•
 | E |  | 4 |   4
 •---•  •---•

 
 Remove C.
 
 'C' is accessible through id 2, if we remove it it will become.
 (swapping both array and indices lut)
 
 array         ids
 •---•  •---•
 | A |  | 0 |   0
 •---•  •---•
 | B |  | 1 |   1
 •---•  •---•
 | E |  | 4 |   2
 •---•  •---•
 | D |  | 3 |   3
 •---•  •---•
 |nul|  | 2 |   4
 •---•  •---•

 
 the id 2 points to the 4'th index, and the data 'E' can still be accessible via
 id 4. even if we move the element at the center, the ids are still mapped to its actual value.
 
 and the id can easily be reused by using free list
*/

zp_cold int zp_container_init(zp_container *const zp_restrict c, const uint16_t stride, const size_t reserve, const float growth_base) {
 assert(reserve != 0);
 assert(stride != 0);
 
 /* not exact, but still it can grow in a decent way, and able to minimize the realloc */
 c->_growth_base = zp_log2(growth_base);
	c->_increase_count = (size_t)zp_ceil(zp_log2((float)growth_base) / c->_growth_base);

	c->_reserve = reserve;
	c->_max_size = reserve;
	c->_stride = stride;
	c->_size = 0;
	
	const size_t bytes_size = c->_max_size * c->_stride;
	const size_t helper_size = sizeof(zp_container_id) * c->_max_size;
	
	const size_t luts_size = helper_size * 3 + bytes_size;
 uint8_t *ptr = (uint8_t*)malloc(luts_size + ZP_MEMORY_ALIGNMENT * 3);
 c->_bytes = ptr;
 ptr += bytes_size;
 c->_to_index_lut = (zp_container_id*)zp_alignto((uint64_t)ptr, ZP_MEMORY_ALIGNMENT);
 ptr = ((uint8_t*)c->_to_index_lut) + helper_size;
 c->_to_id_lut = (zp_container_id*)zp_alignto((uint64_t)ptr, ZP_MEMORY_ALIGNMENT);
 ptr = ((uint8_t*)c->_to_id_lut) + helper_size;
 c->_free_list = (zp_container_id*)zp_alignto((uint64_t)ptr, ZP_MEMORY_ALIGNMENT);

// assert(((uint64_t)(((uint8_t*)c->_free_list) + helper_size) - (uint64_t)c->_bytes) == (luts_size + ZP_MEMORY_ALIGNMENT * 3) );
 assert((((uint64_t)c->_bytes) % 2) == 0);
 assert((((uint64_t)c->_to_index_lut) % ZP_MEMORY_ALIGNMENT) == 0);
 assert((((uint64_t)c->_to_id_lut) % ZP_MEMORY_ALIGNMENT) == 0);
 assert((((uint64_t)c->_free_list) % ZP_MEMORY_ALIGNMENT) == 0);

 if(zp_unlikely(!c->_bytes))
  return -1;
  
 size_t i = c->_size;
 while(i < c->_max_size) {
  c->_to_index_lut[i]._val = i;
  c->_free_list[i]._val = i;
  i++;
 }
 return 0;
}


zp_cold void zp_container_destroy(zp_container *const zp_restrict c) {
 free(c->_bytes);
}

/*
 just return last id, it always points to the last element.

 we can free list to reuse the id that's been freed previously
*/
zp_container_id zp_container_acquire(zp_container *const zp_restrict c) {
 if(c->_size >= c->_max_size)
  if(resize(c))
   return ZP_CONTAINER_NULL_ID;
 zp_container_id allocated = c->_free_list[c->_size];
 c->_to_id_lut[c->_to_index_lut[allocated._val]._val] = allocated;
 c->_size++;
 return allocated;
}

/*
 swap & pop, this creates no gap and as fast as array iteration
 and retrieve the id to the free list
*/
void zp_container_release(zp_container *const zp_restrict c, const zp_container_id id) {
 assert(c->_to_index_lut[id._val]._val < c->_size);
 assert(c->_size != 0);
 assert(!zp_container_id_isnull(id));
 
 size_t current_index = c->_to_index_lut[id._val]._val;
 size_t last_index = c->_size - 1;
 zp_container_id last_id = c->_to_id_lut[last_index];
 
 if(last_index != current_index) {
  memcpy(c->_bytes + ((size_t)current_index * (size_t)c->_stride), c->_bytes + ((size_t)last_index * (size_t)c->_stride), c->_stride);
  memswap(&c->_to_index_lut[id._val], &c->_to_index_lut[last_id._val], sizeof(zp_container_id));
  memswap(&c->_to_id_lut[current_index], &c->_to_id_lut[last_index], sizeof(zp_container_id));
 }
 c->_free_list[last_index] = id;
 c->_size--;
}


zp_cold void zp_container_insertion_sort(zp_container *const zp_restrict c, int (*should_swap)(void*, void*)) {
 for(size_t i = 1; i < c->_size; i++) {
  int j = i;
  int jm1 = j - 1;
  void *j_data = c->_bytes + (j * c->_stride);
  void *jm1_data = c->_bytes + (jm1 * c->_stride);
  
  while (j > 0 && should_swap(jm1_data, j_data)) {
   memswap(j_data, jm1_data, c->_stride);
   zp_container_id j_id = c->_to_id_lut[j];
   zp_container_id jm1_id = c->_to_id_lut[jm1];
   
   memswap(&c->_to_id_lut[j], &c->_to_id_lut[jm1], sizeof(zp_container_id));
   memswap(&c->_to_index_lut[j_id._val], &c->_to_index_lut[jm1_id._val], sizeof(zp_container_id));
   j--;
  }
 }
}


/*FOR DEBUGGING ONLY*/
/*

void is_valid(zp_container* c) {
 for(zp_container_id i = 0; i < c->_size; i++) {
  assert(c->_to_index_lut[c->_to_id_lut[i]] == i);
 }

 for(zp_container_id id = 0; id < c->_size; id++) {
  assert(c->_to_id_lut[c->_to_index_lut[id]] == id);
 }
}

*/



