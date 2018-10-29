/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  11 November 2017                                       */
/*      Purpose:  custom 1D/2D signal processing stuff					 */
/*      Authors:                                                         */
/*        Iurii Efimov                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#ifndef __IMAGEPROCESSING_H__
#define __IMAGEPROCESSING_H__

#include "GlobalParams.h"
#include "BorderDetection.h"

#define fPI         (3.141593f)
#define dPI         (3.14159265358979323846)
#define masekPI     (3.14159265358979)
#ifndef PI
#define PI          fPI
#endif


/************************* HISTOGRAM FUNCTIONS *****************************/
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
	int hw);          // IN:  half-window


void IPL_HIST_Blur_double(
    double* dst,            // OUT: blurred hist
    const double* val,      // IN:  source hist
    int len,                // IN:  length of the hist
    int hw);                // IN:  half-window


int IPL_HIST_mdn_pixE(
	const unsigned char* bufin, // input array
	int                  lstr);  // size bufin

/************************* DYNAMIC TIME WARPING STUFF *****************************/
typedef struct sDTWInstance {
	sPnt path;
	float *mat;
	int sizeL, sizeR;
} sDTWInstance;

RESULT_CODE calcFullDynamicTimeWarpingWrapper(sPnt* pathDTW, int* sizePath, float* matDTW,
	int* pDataL, int sizeL, int* pDataR, int sizeR, int nCompare, float regWeight);

RESULT_CODE myDynamicTimeWarpingPyramid(sPnt* pathDTW, int* sizePath, float* matDTW,
	int* pDataL, int sizeL, int* pDataR, int sizeR, int sizeSeriesWindow, int sizeMetricWindow, int* nScales);

void downsampleSeries(int* dst, int dstSize, int* src, int srcSize, int factor);
void upsamplePath(sPnt* path, int* sizeP, int scaleFactor);

/************************ PROJECTION PROCESSING ***********************************/
RESULT_CODE sequenceGradient(int *pGrData, int *pData, int size);


RESULT_CODE myDynamicTimeWarping(
	sPnt* pathDTW,				// IN/OUT: dynamic time warping path
	float* matDTW,				// IN/OUT: dynamic time warping matrix
	int* sizePath,				// IN/OUT: dynamic time warping path size
	int* pDataL,
	int sizeL,
	int* pDataR,
	int sizeR,
	int sizeSeriesWindow,
	int sizeMetricWindow
);

EXTERN_DLL_EXPORT RESULT_CODE FindHoughProjection(
	int* pProjR,                // OUT: circular projection - right side nums
	int* pProjL,                // OUT: circular projection - left side nums
	int angBeg,                 // IN: start angle for projection
	int angEnd,                 // IN: end angle for projection
	int* Kernel3x3,             // IN: edge detection kernel
	int thrHoriz,               // IN: horizontal derivative threshold
	float thrDot,                 // IN: circular shape threshold
	int minR,
	int maxR,
	unsigned char* im,          // IN: source image
	int xs,                     // IN: image size
	int ys,                     //
	int xc,                     // IN: center x coordinate
	int yc,						// IN: center y coordinate
	int blur_hw);               // blurring half-window size

RESULT_CODE IPL_PROJ_FindHoughDonatorProjection7(
	int* pnProjRN_b,          // OUT: circular projection - right side nums
	int* pnProjLN_b,          // OUT: circular projection - left side nums
	const unsigned char* im,  // IN:  image
	int xs,                   // IN:  image size
	int ys,                   // 
	int xc,                   // IN:  ring center
	int yc,                   //
	int r,                    // IN:  inner radius
	int* pR,                  // IN/OUT: outer radius proposed/actually used
	int nMaxRadOfDence,       // IN:  maximum radius for dense processing. 0 - dense always
	void* buf,                // IN:  external buffer
	int* buflen);              // IN/OUT: allocated/used bytes

RESULT_CODE IPL_PROJ_FindHoughDonatorProjection6(
	int* pnProjRN_b,          // OUT: circular projection - right side nums
	int* pnProjLN_b,          // OUT: circular projection - left side nums
	const unsigned char* im,  // IN:  image
	int xs,                   // IN:  image size
	int ys,                   // 
	int xc,                   // IN:  ring center
	int yc,                   //
	int r,                    // IN:  inner radius
	int* pR,                  // IN/OUT: outer radius proposed/actually used
	int nMaxRadOfDence,       // IN:  maximum radius for dense processing. 0 - dense always
	void* buf,                // IN:  external buffer
	int* buflen);              // IN/OUT: allocated/used bytes

RESULT_CODE IPL_PROJ_FindHoughDonatorProjection5(
	int* pnProjRN_b,          // OUT: circular projection - right side nums
	int* pnProjLN_b,          // OUT: circular projection - left side nums
	const unsigned char* im,  // IN:  image
	int xs,                   // IN:  image size
	int ys,                   // 
	int xc,                   // IN:  ring center
	int yc,                   //
	int r,                    // IN:  inner radius
	int* pR,                  // IN/OUT: outer radius proposed/actually used
	int nMaxRadOfDence,       // IN:  maximum radius for dense processing. 0 - dense always
	void* buf,                // IN:  external buffer
	int* buflen);

// == image filtering functions ==============================================

// 3x3 Haussian blur (1-2-1)
EXTERNC void IPL_FILT_HaussBlur3x3(
	unsigned char* dst,         // OUT: dst image
	const unsigned char* src,   // IN:  src image
	int xs,                     // IN:  image size
	int ys);                    //

// integer 3x3 Haussian blur (1-2-1)
EXTERNC void IPL_FILT_HaussBlur3x3_Int(
	int* dst,                   // OUT: dst image
	const int* src,             // IN:  src image
	int xs,                     // IN:  image size
	int ys,                     //
	int mode);  // 0 - normalized, 1- no normalization

				// Blur an image with a gaussian filter
EXTERNC RESULT_CODE IPL_FILT_GaussianSmooth_ext_uint8(
	short *smoothedim,
	double *tempim,
	const unsigned char *image,
	int xs,
	int ys,
	double sigma_x,
	double sigma_y,
	double _BOOSTBLURFACTOR);


void Dilate3x3Cross(
	uint8* dst,
	const uint8* src,
	int W,
	int H);

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

/******************** EDGE DETECTION ******************/
// canny edge detection
EXTERNC RESULT_CODE IPL_FILT_Canny(
	unsigned char *edge,        // OUT: canny mask
	const unsigned char *image, // IN:  source image
	short ** ppgx,              // OUT: map of horizontal gradients
	short ** ppgy,              // OUT: map of vertical gradients
	int xs,                     // IN:  image size
	int ys,                     // 
	double sigma,               // IN:  haussian sigma
	double tlow,                // IN:  low hysteresis threshold (turning off)
	double thigh,               // IN:  high hysteresis threshold (turning on)
	void* pvBuf,                // IN:  temporary buffer
	int* pnLen,                 // IN/OUT: bytes allocated/used
	const char *nambeg);        // IN:  set this to NULL
#endif