#include "zp_physics/core2d/contacthash2d.h"
#include "zp_physics/math.h"
#include <stdlib.h>
#include <string.h>

/*
Robert Jenkins' 32 bit integer hash function
 src : https://gist.github.com/badboy/6267743
*/
zp_const zp_inline uint32_t hash(uint32_t a) {
 a = (a + 0x7ed55d16) + (a << 12);
 a = (a ^ 0xc761c23c) ^ (a >> 19);
 a = (a + 0x165667b1) + (a << 5);
 a = (a + 0xd3a2646c) ^ (a << 9);
 a = (a + 0xfd7046c5) + (a << 3);
 a = (a ^ 0xb55a4f09) ^ (a >> 16);
 return a;
}

zp_const zp_inline zp_hot uint32_t make_hash_id(uint16_t a, uint16_t b) {
 if(a > b) {
  uint16_t tmp = b;
  b = a;
  a = tmp;
 }
 uint32_t combined_id = ((((uint32_t)a) << 16) | (uint32_t)b);
 return hash(combined_id);
}

/*
 same as insert function, but no rehash to avoid recursion.
*/
zp_noinline zp_hot static void rehash_insert(zp_contacthash2d *const zp_restrict hash, const zp_manifold2d *const zp_restrict m) {
 const size_t hash_index = make_hash_id(m->_body_a._val, m->_body_b._val) & hash->_bucket_index_mask;

 zp_contacthash2d_node *start_node = NULL;

 zp_container_id next_index = hash->_bucket[hash_index];
 size_t it = 0;


 while(!zp_container_id_isnull(next_index)) {
  assert((it < hash->_memory_pool._size) && "infinite loop.");
  assert((hash->_memory_pool._to_index_lut[next_index._val]._val < hash->_memory_pool._size) && "invalid id");
 	start_node = (zp_contacthash2d_node*)zp_container_get(&hash->_memory_pool, next_index);

 	if((zp_container_id_isequal(start_node->_value._body_a, m->_body_a) && zp_container_id_isequal(start_node->_value._body_b, m->_body_b)) || (zp_container_id_isequal(start_node->_value._body_a, m->_body_b) && zp_container_id_isequal(start_node->_value._body_b, m->_body_a))) {
 	 start_node->_queried = 1;
 	 zp_manifold2d_combine(&start_node->_value, m);
 	 return;
 	}
 	it++;
		next_index = start_node->_next;
 }
  
	next_index = zp_container_acquire(&hash->_memory_pool);
 assert((hash->_memory_pool._to_index_lut[next_index._val]._val < hash->_memory_pool._size) && "invalid id");
	start_node = (zp_contacthash2d_node*)zp_container_get(&hash->_memory_pool, next_index);
	start_node->_queried = 2;
	start_node->_next = hash->_bucket[hash_index];
	start_node->_prev = ZP_CONTAINER_NULL_ID;
	start_node->_allocation = next_index;
 memcpy(&start_node->_value, m, sizeof(zp_manifold2d));
 
 if(!zp_container_id_isnull(hash->_bucket[hash_index])) {
 	assert((hash->_memory_pool._to_index_lut[hash->_bucket[hash_index]._val]._val < hash->_memory_pool._size) && "invalid id");
 	start_node = (zp_contacthash2d_node*)zp_container_get(&hash->_memory_pool, hash->_bucket[hash_index]);
	 start_node->_prev = next_index;
 }
 hash->_bucket[hash_index] = next_index;
}



/*
 separate chaining combined with dense container Hash Set
 this able to handle large data loads, no cluster issue and 
 somehow cache friendly 
*/
zp_cold int zp_contacthash2d_init(zp_contacthash2d *const zp_restrict hash, const size_t bucket_size, const size_t reserve, const float growth_base) {
 hash->_bucket_size = zp_nextp2(bucket_size);
 hash->_bucket_index_mask = hash->_bucket_size - 1;
 hash->_bucket = (zp_container_id*)malloc(sizeof(zp_container_id) * hash->_bucket_size);
 memset(hash->_bucket, 0xFF, sizeof(zp_container_id) * hash->_bucket_size);

 if(zp_unlikely(zp_container_init(&hash->_memory_pool, sizeof(zp_contacthash2d_node), reserve, growth_base)))
  return -1;
 return 0;
}


zp_cold void zp_contacthash2d_destroy(zp_contacthash2d *const zp_restrict hash) {
 free(hash->_bucket);
 zp_container_destroy(&hash->_memory_pool);
}


void zp_contacthash2d_insert(zp_contacthash2d *const zp_restrict hash, const zp_manifold2d *const zp_restrict m) {
 if(((float)hash->_memory_pool._size / (float)hash->_bucket_size) > 0.5f) {
  /* rehash, this maintains O(1) average */
  zp_contacthash2d new_hash;
  if(zp_unlikely(zp_contacthash2d_init(&new_hash, hash->_bucket_size * 2, hash->_memory_pool._reserve * 2, hash->_memory_pool._growth_base)) ) {
   return;
  }
  for(uint16_t i = 0; i < hash->_memory_pool._size; i++) {
   zp_contacthash2d_node *current_node = (zp_contacthash2d_node*)(hash->_memory_pool._bytes + ((size_t)i * (size_t)hash->_memory_pool._stride));
   rehash_insert(&new_hash, &current_node->_value);
  }
  zp_contacthash2d_destroy(hash);
  memcpy(hash, &new_hash, sizeof(zp_contacthash2d));
 }
 
 const size_t hash_index = make_hash_id(m->_body_a._val, m->_body_b._val) & hash->_bucket_index_mask;

 zp_contacthash2d_node *start_node = NULL;

 zp_container_id next_index = hash->_bucket[hash_index];
 size_t it = 0;


 while(!zp_container_id_isnull(next_index)) {
  assert((it < hash->_memory_pool._size) && "infinite loop.");
  assert((hash->_memory_pool._to_index_lut[next_index._val]._val < hash->_memory_pool._size) && "invalid id");
 	start_node = (zp_contacthash2d_node*)zp_container_get(&hash->_memory_pool, next_index);

 	if((zp_container_id_isequal(start_node->_value._body_a, m->_body_a) && zp_container_id_isequal(start_node->_value._body_b, m->_body_b)) || (zp_container_id_isequal(start_node->_value._body_a, m->_body_b) && zp_container_id_isequal(start_node->_value._body_b, m->_body_a))) {
 	 start_node->_queried = 1;
 	 zp_manifold2d_combine(&start_node->_value, m);
 	 return;
 	}
 	it++;
		next_index = start_node->_next;
 }
 
	next_index = zp_container_acquire(&hash->_memory_pool);
 assert((hash->_memory_pool._to_index_lut[next_index._val]._val < hash->_memory_pool._size) && "invalid id");
	start_node = (zp_contacthash2d_node*)zp_container_get(&hash->_memory_pool, next_index);
	start_node->_queried = 1;
	start_node->_next = hash->_bucket[hash_index];
	start_node->_prev = ZP_CONTAINER_NULL_ID;
	start_node->_allocation = next_index;
 memcpy(&start_node->_value, m, sizeof(zp_manifold2d));
 
 if(!zp_container_id_isnull(hash->_bucket[hash_index])) {
 	assert((hash->_memory_pool._to_index_lut[hash->_bucket[hash_index]._val]._val < hash->_memory_pool._size) && "invalid id");
 	start_node = (zp_contacthash2d_node*)zp_container_get(&hash->_memory_pool, hash->_bucket[hash_index]);
	 start_node->_prev = next_index;
 }
 hash->_bucket[hash_index] = next_index;
}



/*
 if its not queried, more likely collision does not exist anymore,
 so remove them.
*/
void zp_contacthash2d_remove_unused(zp_contacthash2d *const zp_restrict hash) {
 uint16_t i = 0;
 while(i < hash->_memory_pool._size) {
 
  zp_contacthash2d_node *current_node = (zp_contacthash2d_node*)(hash->_memory_pool._bytes + ((size_t)i * (size_t)hash->_memory_pool._stride));

  if(current_node->_queried) {
   current_node->_queried = 0;
  } else {
   if(!zp_container_id_isnull(current_node->_prev)) {
     assert((hash->_memory_pool._to_index_lut[current_node->_prev._val]._val < hash->_memory_pool._size) && "invalid id");
    ((zp_contacthash2d_node*)zp_container_get(&hash->_memory_pool, current_node->_prev))->_next = current_node->_next;
   } else
    hash->_bucket[make_hash_id(current_node->_value._body_a._val, current_node->_value._body_b._val) & hash->_bucket_index_mask] = current_node->_next;
  
   if(!zp_container_id_isnull(current_node->_next)) {
    assert((hash->_memory_pool._to_index_lut[current_node->_next._val]._val < hash->_memory_pool._size) && "invalid id");
    ((zp_contacthash2d_node*)zp_container_get(&hash->_memory_pool, current_node->_next))->_prev = current_node->_prev;
   }
    
   assert((hash->_memory_pool._to_index_lut[current_node->_allocation._val]._val < hash->_memory_pool._size) && "invalid id");
   zp_container_release(&hash->_memory_pool, current_node->_allocation);
   continue;
  }
  i++;
 }
}


/*
 TODO : needs to be fast, since its used to check if its recently queried 
 and eliminate the duplicates.
*/
uint8_t zp_contacthash2d_is_queried(zp_contacthash2d *const zp_restrict hash, const zp_container_id a, const zp_container_id b) {
 const size_t hash_index = make_hash_id(a._val, b._val) & hash->_bucket_index_mask;
 zp_contacthash2d_node *start_node = NULL;

 zp_container_id next_index = hash->_bucket[hash_index];
 size_t it = 0;


 while(!zp_container_id_isnull(next_index)) {
  assert((it < hash->_memory_pool._size) && "infinite loop.");
  assert((hash->_memory_pool._to_index_lut[next_index._val]._val < hash->_memory_pool._size) && "invalid id");
 	start_node = (zp_contacthash2d_node*)zp_container_get(&hash->_memory_pool, next_index);

 	if((zp_container_id_isequal(start_node->_value._body_a, a) && zp_container_id_isequal(start_node->_value._body_b, b)) || (zp_container_id_isequal(start_node->_value._body_a, b) && zp_container_id_isequal(start_node->_value._body_b, a))) {
	  return start_node->_queried;
 	}
 	it++;
		next_index = start_node->_next;
 }
 return 0xFF;
}

