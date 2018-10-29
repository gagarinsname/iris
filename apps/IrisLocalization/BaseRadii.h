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

/*RESULT_CODE GenerateProjectionFeatures(
    char* srcFolder,            // source images folder
    char* srcMarkupFile,        // markup file for them
    char* dstFile);             // resulting feature file
*/


RESULT_CODE myDetectBaseRadii(
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

int DetectBaseRadii(
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
