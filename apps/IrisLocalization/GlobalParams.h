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
#include "stddefs.h"

#define EXTERN_DLL_EXPORT extern __declspec(dllexport)  

#define RESULT_CODE    int
enum ERROR_CODES {
    ERROR_OK,
    ERROR_MEMORY,
    ERROR_NULL_POINTER,
    ERROR_WRONG_INPUT,
	ERROR_NO_ANSWER
};

#define TIMERS_ENABLED					(0)
#define INTSQRT(x) (int)(sqrtf((float)x))

#define int40 long long
typedef unsigned char uint8;

#define MAX_FILENAME 1024


#define MIA_abs(x) (((x)>0)?(x):(-(x)))
#define MIA_bound(x,a,b) (((x)<(a))?(a):(((x)>(b))?(b):(x)))
#define MIA_roundf(x) ((int32)(((x)>=0)?((x)+0.5f):((x)-0.5f)));
#define MIA_max(a,b) (((a) > (b)) ? (a) : (b))
#define MIA_min(a,b) (((a) < (b)) ? (a) : (b))

typedef struct
{
	int xc;   // X coordinate
	int yc;   // Y coordinate
	int r;    // radius (elliptic big axis)
} CInfo;

/******** Paired Gradient-based segmentation related stuff **********************************/
#define BDL_PGRAD_CALCRADIUS    0x00000001  // perform radius calculation
#define BDL_PGRAD_SAVEACC		0x00000002	// save the accumulator image
#define BDL_PGRAD_SAVECANNY		0x00000004	// save the canny edge image
#define BDL_PGRAD_LINEVOTE		0x00000008	// voting in the accumulator using lines
#define BDL_PGRAD_USE_CANNY		0x00000010	// use canny to detect edges
#define BDL_PGRAD_USE_SOBEL		0x00000020	// use sobel to detect edges
#define BDL_PGRAD_USE_IMAGE		0x00000040	// use edgemap from input
#define BDL_PGRAD_CANNY_AUTO	0x00000100	// auto selection for canny params

#define BDL_PUPREF_SAVECROP		0x00000020	//save the cropped border image for the pupil
#define BDL_PUPREF_SAVEPOLAR	0x00000040	//save the polar transformed pupil border region
#define FLT_SOBEL_SAVEIMAGE		0x00000080	//save the result of filter applying as an image
#define FLT_SOBEL_SELECTION		0x00000100	//select pixels with the top 15% gradient magnitudes
#define BDL_PGRAD_PREPROCESSING 0x00000200

#define MINIMAL_PUPIL_RADIUS 50
#define BLUR_HIST_WND 4

#endif