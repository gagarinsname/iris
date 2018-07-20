/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  11 November 2017                                       */
/*      Purpose:  global parameters for base radius detection API        */
/*      Authors:                                                         */
/*        Iurii Efimov                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#ifndef __GlobalParams_h__
#define __GlobalParams_h__

#include <math.h>

#define RESULT_CODE    int

enum ERROR_CODES {
    ERROR_OK,
    ERROR_MEMORY,
    ERROR_NULL_POINTER,
    ERROR_WRONG_INPUT
};

#define INTSQRT(x) (int)(sqrtf((float)x))

#define int40 long long
typedef unsigned char uint8;

#define MAX_FILENAME 1024

#define MIA_abs(x) (((x)>0)?(x):(-(x)))
#define MIA_bound(x,a,b) (((x)<(a))?(a):(((x)>(b))?(b):(x)))
#define MIA_roundf(x) ((int32)(((x)>=0)?((x)+0.5f):((x)-0.5f)));
#define MIA_max(a,b) (((a) > (b)) ? (a) : (b))
#define MIA_min(a,b) (((a) < (b)) ? (a) : (b))

#endif