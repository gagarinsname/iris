/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  20 March 2001                                          */
/*      Revision: 1.0.00                                                 */
/*      Purpose:  Standard defines used everywhere in MIA soft           */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*        Alexander Gneushev (added)                                                               */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#ifndef __stddefs_h__
#define __stddefs_h__

#include <memory.h>
#include "unitypes.h"

#ifdef _MSC_VER
  #pragma warning(disable: 4100)  // unreferenced formal parameter
#endif

/* frequently used: declaration for exporting to c++ */
#ifndef EXTERNC
  #ifdef __cplusplus
    #define EXTERNC extern "C"
  #else
    #define EXTERNC
  #endif
#endif

/* frequently used: null */
#if !defined(NULL)
  #if defined(__cplusplus)
    #define NULL 0L
  #else
    #define NULL ((void*)0)
  #endif
#endif

/* frequently used: unique type of the handle */
#ifndef DECL_HANDLE
  #define DECL_HANDLE(name) typedef struct { uint32 unused; } *name
#endif

/* frequently used: making four-character IDs for data blocks */
#ifndef GENERATE_FOURCC_ID
  #define GENERATE_FOURCC_ID(ch0, ch1, ch2, ch3) \
              ( ((uint32)(uint8)(ch0) <<  0) | \
                ((uint32)(uint8)(ch1) <<  8) | \
                ((uint32)(uint8)(ch2) << 16) | \
                ((uint32)(uint8)(ch3) << 24)   )
#endif

/* Static declaration for functions */
#ifdef staticFunc
  #error staticFunc should be defined once in this file
#endif
#define staticFunc
/*#define staticFunc static */

/* Static declaration for variables */
#ifdef staticVar
  #error staticVar should be defined once in this file
#endif
#define staticVar
/*#define staticVar static */

/* Static declaration for internal function variables */
#ifdef staticInt
  #error staticInt should be defined once in this file
#endif
#define staticInt static

/* Static declaration for constants */
#ifdef staticConst
  #error staticConst should be defined once in this file
#endif
/*#define staticConst const*/
#define staticConst static const

/* memory management functions */
#define MIA_memset memset
#define MIA_memcpy memcpy

/* math functions */
#define MIA_abs(x) (((x)>0)?(x):(-(x)))
#define MIA_bound(x,a,b) (((x)<(a))?(a):(((x)>(b))?(b):(x)))
#define MIA_roundf(x) ((int32)(((x)>=0)?((x)+0.5f):((x)-0.5f)));
#define MIA_max(a,b) (((a) > (b)) ? (a) : (b))
#define MIA_min(a,b) (((a) < (b)) ? (a) : (b))

#define fPI         (3.141593f)
#define dPI         (3.14159265358979323846)
#define masekPI     (3.14159265358979)
#ifndef PI
  #define PI          fPI
#endif

#endif /*__stddefs_h__*/
