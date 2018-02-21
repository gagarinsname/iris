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

#define INTSQRT(x) (int)(sqrt(x))

#define int40 long long
typedef unsigned char uint8;

#define MAX_FILENAME 1024

#endif