/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  11 November 2017                                       */
/*      Purpose:  internal base radius detection API                     */
/*      Authors:                                                         */
/*        Ivan Matveev, Iurii Efimov                                     */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#ifndef __BaseRadii_h__
#define __BaseRadii_h__

#include "GlobalParams.h"
#include "BorderDetection.h"

#define MINIMAL_PUPIL_RADIUS 50
#define EXTERN_DLL_EXPORT extern __declspec(dllexport)  

typedef struct SMarkup
{
    char filename[MAX_FILENAME];
    SCenterInfo pupil, iris;
    SBazRadInfo BazRad;
} SMarkup;


/*RESULT_CODE GenerateProjectionFeatures(
    char* srcFolder,            // source images folder
    char* srcMarkupFile,        // markup file for them
    char* dstFile);             // resulting feature file
*/

/************************* DYNAMIC TIME WARPING STUFF *****************************/
typedef struct sDTWInstance {
	sPnt path;
	float *mat;
	int sizeL, sizeR;
} sDTWInstance;

EXTERN_DLL_EXPORT RESULT_CODE calcFullDynamicTimeWarpingWrapper(sPnt* pathDTW, int* sizePath, float* matDTW,
	int* pDataL, int sizeL, int* pDataR, int sizeR, int nCompare, float regWeight);

EXTERN_DLL_EXPORT RESULT_CODE myDynamicTimeWarpingPyramid(sPnt* pathDTW, int* sizePath, float* matDTW,
	int* pDataL, int sizeL, int* pDataR, int sizeR, int sizeSeriesWindow, int sizeMetricWindow, int* nScales);

EXTERN_DLL_EXPORT void downsampleSeries(int* dst, int dstSize, int* src, int srcSize, int factor);
EXTERN_DLL_EXPORT void upsamplePath(sPnt* path, int* sizeP, int scaleFactor);

/************************ PROJECTION PROCESSING ***********************************/
EXTERN_DLL_EXPORT RESULT_CODE sequenceGradient(int *pGrData, int *pData, int size);


EXTERN_DLL_EXPORT void IPL_HIST_Blur(
	int* dst,         // OUT: blurred hist
	const int* val,   // IN:  source hist
	int len,          // IN:  length of the hist
	int hw);           // IN:  half-window

EXTERN_DLL_EXPORT RESULT_CODE myDynamicTimeWarping(
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

EXTERN_DLL_EXPORT RESULT_CODE myDetectBaseRadii(
    SPupilInfo* psPI,         // OUT: pupil data
    SIrisInfo* psII,          // OUT: iris data
    SCenterInfo* psCI,        // IN:  center data
    unsigned char* im,        // IN:  source image
    int W,                    // IN:  image size
    int H,                    //
    int* kernel3x3,           // edge detection kernel: kernel3x3[] = { -3, 0, 3, -10, 0, 10, -3, 0, 3 };
    int angMin,               // minimal projection estimation angle
    int angMax,               // maximal projection estimation angle
    int minR,                 // minimal radius
    int maxR,                 // maximal radius
    int thrHoriz,             // threshold for horizontal derivative
    float thrDot,             // threshold for circular shape: 0.8 - 1.0
    int* pProjR,              // circular projection - right side nums
    int* pProjL,              // circular projection - left side nums
	int blurHW
);

#define BLUR_HIST_WND 4
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

EXTERN_DLL_EXPORT int DetectBaseRadii(
    SPupilInfo* psPI,         // OUT: pupil data
    SIrisInfo* psII,          // OUT: iris data
    const SCenterInfo* psCI,  // IN:  center data
    const unsigned char* im,  // IN:  source image
    int xs,                   // IN:  image size
    int ys,                   //
    int* pProjR,               // circular projection - right side nums
    int* pProjL               // circular projection - left side nums
    // int mode // mask: 1- draw projs to buf, 2- use projs from buf
);

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

// detect approximate pupil and iris by projection method
RESULT_CODE IVIR_PrPu(
    SPupilInfo* psPI,         // OUT: pupil data
    SIrisInfo* psII,          // OUT: iris data
    const SCenterInfo* psCI,  // IN:  center data
    const unsigned char* im,  // IN:  source image
    int xs,                   // IN:  image size
    int ys,                   //
    int mode, // mask: 1- draw projs to buf, 2- use projs from buf
    int br_size,    // predefined iris size (-1 - not predefined)
    int br_spread,  // predefined spread of iris size (-1 - not predefined)
    int* projR,
    int* projL,
    const char* nambeg);       // used for debug, set to NULL

// detection of bazrad with red-eye effect
RESULT_CODE IVIR_PrPuRE(
	SBazRadInfo* psBI,        // OUT: bazrad data
	const SCenterInfo* psCI,  // IN:  center data
	const unsigned char* im,  // IN:  source image
	int xs,                   // IN:  image size
	int ys,                   //
	int mode, // mask: 1- draw projs to buf, 2- use projs from buf
	void* buf,                // IN:  temporary buffer
	int* pnLoc,               // OUT: buffer size outputted
	const char* nambeg);      // used for debug, set to NULL

// Estimates iris X coord and radius. Sets Y equal to pupil. 
// Outputs two qualities: q1 for right side and q2 for left side. 
// pIris->OccUpB is set to special flag
// pIris->OccUpE is set to maximum iris radius searched
RESULT_CODE IVIR_FindApproximateIrisByPupil_br(
	SIrisInfo* pIris,           // OUT: iris data
	const SPupilInfo* pPupil,   // IN:  pupil data
	const unsigned char* im,    // IN:  source image
	int xs,                     // IN:  image size
	int ys,                     //
	int nMaxRadOfDence,         // IN:  maximum radius for dence processing. 0 - dence always
	void* buf,                  // IN:  buffer
	int* buflen,                // IN/OUT: allocated/used buffer size
	const char* nambeg);        // IN:  used for debugging, set it to NULL

								// Estimates iris X coord and radius. Sets Y equal to pupil. 
RESULT_CODE IVIR_FindApproximateIrisByPupil2(
	SIrisInfo* pIris_,          // OUT: iris data
	const SPupilInfo* pPupil,   // IN:  pupil data
	const unsigned char* im,    // IN:  source image
	int xs,                     // IN:  image size
	int ys,                     //
	int nMaxRadOfDence,         // IN:  maximum radius for dence processing. 0 - dence always
	void* buf,                  // IN:  buffer
	int* buflen,                // IN/OUT: allocated/used buffer size
	const char* nambeg);        // IN:  used for debugging, set it to NULL

								// Estimates iris X coord and radius. Sets Y equal to pupil. 
								// Outputs two qualities: q1 for right side and q2 for left side. 

RESULT_CODE IVIR_FindApproximateIrisByPupil3(
	SIrisInfo* pIris,           // OUT: iris data
	const SPupilInfo* pPupil,   // IN:  pupil data
	const unsigned char* im,    // IN:  source image
	int xs,                     // IN:  image size
	int ys,                     //
								//  int nFlags,                 // IN:  processing mode. 1-generate projections, 2 - use projs.
	int* pnProjs,               // IN/OUT: projections
	void* buf,                  // IN:  temporary buffer
	int* buflen,                // IN/OUT: allocated/used buffer size
	const char* nambeg);

// same as '3' but with awareness of left/right eye
RESULT_CODE IVIR_FindApproximateIrisByPupil4(
	SIrisInfo* pIris,           // OUT: iris data
	const SPupilInfo* pPupil,   // IN:  pupil data
	const unsigned char* im,    // IN:  source image
	int xs,                     // IN:  image size
	int ys,                     //
	int nFlags,                 // IN:
	int* pnProjs,               // IN/OUT: projections
	void* buf,                  // IN:  temporary buffer
	int* buflen,                // IN/OUT: allocated/used buffer size
	const char* nambeg);

// based on pupil info, find projection local maxima positions suitable for iris
RESULT_CODE IVIR_FindMaximaPositions(
	int* pMaxes,        // OUT: positions of local maxima (sorted by value)
	int* vMaxes,        // OUT: values of local maxima
	int** pProjL,       // OUT: pointer to left projection buffer (out of pvbuf)
	int** pProjR,       // OUT: pointer to right projection buffer (out of pvbuf)
	int* pnProjLNorm,   // OUT: left projection norm
	int* pnProjRNorm,   // OUT: right projection norm
	int r,              // IN:  inner radius
	int* pR,            // IN/OUT: outer radius
	const unsigned char* im,  // IN:  source image
	int xs,                   // IN:  image size
	int ys,                   //
	const SPupilInfo* pPupil, // IN:  pupil data
	int nMaxes,               // IN:  number of required maxima
	void* pvBuf,              // IN:  temporary buffer
	int* pnLen,               // IN/OUT: number of bytes allocated/used
	const char* nambeg);      // IN:  for debugging purposes, set to NULL

int IVIR_ProjPostProcess(  // returns unnormed histogram maximum
	int* hist,
	int len);

int ownIVIR_ListMaxs(
	int* nPos,
	const int* hist1,
	int minrad,
	int len);

int IVIR_P2I_GeomQuality(
	int PupilR,
	int posL,
	int posR);

#endif //__BaseRadii_h__
