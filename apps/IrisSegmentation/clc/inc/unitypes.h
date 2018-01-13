/************************************************************* 
  unitypes.h

  Description: Definition of portable types.
 
  ANG      initial C version (16-bit and 32-bit Windows)
  ANG      added DSP TMS320C6x relevant
  MIA      cut dummies
  MIA      20100528 get rid of 16-bit junk, added 64-bit relevant
 *************************************************************/

#ifndef __unitypes_h__
#define __unitypes_h__

// compile in strict mode
#ifndef STRICT
  #define STRICT
#endif

// 8-bit integers
#ifndef int8
  typedef signed char int8;
#endif

#ifndef uint8
  typedef unsigned char uint8;
#endif

// 16-bit integers
#ifndef int16
  typedef signed short int16;
#endif

#ifndef uint16
  typedef unsigned short uint16;
#endif

// 32-bit integers. Forget about 16-bit platform! 'int' is 32-bit forever!
#ifndef int32
  typedef signed int int32;
#endif

#ifndef uint32
  typedef unsigned int uint32;
#endif

// 64-bit integers
#ifndef int64
  #if (defined(_MSC_VER)||defined(__INTEL_COMPILER))
    typedef signed __int64 int64;
  #else
    // GNU
    typedef signed long long int64;
  #endif
#endif

#ifndef uint64
  #if (defined(_MSC_VER)||defined(__INTEL_COMPILER))
    typedef unsigned __int64 uint64;
  #else
    // GNU
    typedef unsigned long long uint64;
  #endif
#endif

// TMS320C6x does not have 64-bit integer. 
// Instead, for values bigger than 2^32 it has 40-bit and calls it 'long'.
#ifndef int40
  #if defined(_TMS320C6X)
    typedef signed long int40;
  #else
    // non-TMS use 64-bit here
    typedef int64 int40;
  #endif
#endif

#ifndef uint40
  #if defined(_TMS320C6X)
    typedef unsigned long uint40;
  #else
    // non-TMS use 64-bit here
    typedef uint64 uint40;
  #endif
#endif

// type of integer, that has same size as pointer, so as: sizeof(ptr_t)=sizeof(void*)
#ifndef ptr_t
  #if ((defined(__x86_64__))||(defined(_WIN64))||(defined(WIN64)))
    typedef uint64 ptr_t;
  #else
    typedef uint32 ptr_t;
  #endif
#endif

// Fastest possible int type. 16 bits only warrantied.
#ifdef _TMS320C6X
  typedef int16 FASTINT16;
#else
  typedef int32 FASTINT16;
#endif

#endif //__unitypes_h__
