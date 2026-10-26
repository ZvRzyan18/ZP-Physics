#ifndef OXB_MEMORY_H
#define OXB_MEMORY_H

#include "oxb1tz/types.h"

void* oxb_aligned_alloc(oxb_allocator *oxb_restrict const allocator, const size_t bytes_size);
void oxb_aligned_free(oxb_allocator *oxb_restrict const allocator, void *oxb_restrict const ptr);

void oxb_memory_copy(void *oxb_restrict dst, const void *oxb_restrict src, size_t n);
void oxb_memory_set(void *oxb_restrict dst, const uint8_t src, size_t n);

#endif

