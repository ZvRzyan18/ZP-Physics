#ifndef ZP_HINT_H
#define ZP_HINT_H

#include <stddef.h>
#include <stdint.h>


#ifdef __cplusplus
#define zp_restrict __restrict
#else
#define zp_restrict restrict
#endif


#define ZP_MEMORY_ALIGNMENT 16


#if defined(__clang__) || defined(__GNUC__)
#define zp_has_complex_extension

#define zp_likely(x)   __builtin_expect(!!(x), 1)
#define zp_unlikely(x) __builtin_expect(!!(x), 0)
#define zp_inline __attribute__((always_inline, hot)) static inline 
#define zp_noinline __attribute__((noinline))
#define zp_unreachable() __builtin_unreachable()
#define zp_const __attribute__((const))
#define zp_pure __attribute__((pure))
#define zp_hot __attribute__((hot))
#define zp_cold __attribute__((cold))

#define zp_compiler_memory_barrier() __asm__ volatile("" ::: "memory");

#elif defined(_MSC_VER)
#include <intrin.h>


#define zp_likely(x) (x)
#define zp_unlikely(x) (x)
#define zp_inline __forceinline static inline 
#define zp_noinline __declspec(noinline)
#define zp_unreachable() __assume(0)
#define zp_const 
#define zp_pure
#define zp_hot 
#define zp_cold 

#define zp_compiler_memory_barrier() _ReadWriteBarrier()

#else

#define zp_likely(x)
#define zp_unlikely(x)
#define zp_inline
#define zp_noinline
#define zp_unreachable() 
#define zp_const 
#define zp_pure
#define zp_hot 
#define zp_cold 

#define zp_compiler_memory_barrier() 

#endif



//#define ZP_SIMD

/* TODO : implement simd */
#if defined(__ARM_NEON__) || defined(__ARM_NEON) && defined(ZP_SIMD)
#include <arm_neon.h>
#define ZP_NEON

#if !defined(__aarch64__) && !defined(__ARM_FEATURE_FMA)
#define __zp_neon_fma vmla_f32
#define __zp_neon_fmaq vmlaq_f32
#else
#define __zp_neon_fma vfma_f32
#define __zp_neon_fmaq vfmaq_f32
#endif

#endif





#endif

