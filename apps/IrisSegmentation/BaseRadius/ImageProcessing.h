
#ifndef __IMAGEPROCESSING_H__
#define __IMAGEPROCESSING_H__

#include "GlobalParams.h"

#define fPI         (3.141593f)
#define dPI         (3.14159265358979323846)
#define masekPI     (3.14159265358979)
#ifndef PI
#define PI          fPI
#endif


//== histogram functions ====================================================
RESULT_CODE IPL_HIST_CalcMedianInLine(
    unsigned char* line,
    const unsigned char* im,
    int xs,
    int ys,
    int h,
    int kew);
// Returns: median

void IPL_HIST_Blur(
    int* dst,         // OUT: blurred hist
    const int* val,   // IN:  source hist
    int len,          // IN:  length of the hist
    int hw);           // IN:  half-window

void IPL_HIST_Blur_double(
    double* dst,            // OUT: blurred hist
    const double* val,      // IN:  source hist
    int len,                // IN:  length of the hist
    int hw);                // IN:  half-window


// == image filtering functions ==============================================

RESULT_CODE IPL_FILT_GaussianSmooth_uint8(
    short *smoothedim,
    double *tempim,
    const unsigned char *image,
    int xs,
    int ys,
    double sigma);

RESULT_CODE IPL_FILT_GaussianSmooth_ext_uint8(
    short *smoothedim,
    double *tempim,
    const unsigned char *image,
    int xs,
    int ys,
    double sigma_x,
    double sigma_y,
    double _BOOSTBLURFACTOR);

RESULT_CODE IPL_FILT_GaussianSmooth_double(
    double *smoothedim,
    double *tempim,
    const double *image,
    int xs,
    int ys,
    double sigma);

#endif