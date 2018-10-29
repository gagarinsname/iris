/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  16 Feb 2017                                            */
/*      Modified: many times                                             */
/*      Purpose:  Iris border detection functions						 */
/*      Authors:                                                         */
/*        Yuriy Efimov                                                   */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/
#include "stddefs.h"
#include "unitypes.h"
#include "BorderDetection.h"

#define CONTOUR_QUANTIZATION 360
#define BORDER_VALUE 10
#define NORMALIZED_IRIS_HEIGHT 64

#define VOTING_LINES (0)
#define ANGDELT 8
#define RADDELT 5
#define MINRAD 20

/* PUPIL CSP */
#define REFINE_PUPIL 0
#define CDER 50
#define MAXCOST 214748364

typedef struct SEdgeDetectorParams
{
	double tHigh, tLow, tSigma;
} SEdgeDetectorParams;

EXTERNC int IS_Initialize(SSegmentationResult *sResult, char* name);
EXTERNC int IrisSegmentation(SSegmentationResult* Result, uint8 *imgInput, int H, int W, int angle, int flags);
EXTERNC int IS_Deinitialize(SSegmentationResult *sResult);
EXTERNC int IS_ApproxBothBorders(SSegmentationResult *Result, const uint8 *imgInput, int H, int W, int dAngle, int flags); 

EXTERNC RESULT_CODE DetectCenterPairedGradients(int* cenX, int* cenY, SEdgePnt** ppnEdgePnts, int* nGp,
	const uint8* img, int H, int W, float tHigh, float tLow, float sigma, int dPhi, int flags);
EXTERNC RESULT_CODE DetectCenterFromEdgeMapPairedGradients(int* cenX, int* cenY, SEdgePnt** ppnEdgePnts, int* nGp,
	const uint8* img, int H, int W, uint8* pucEdge, int16* pGx, int16* pGy, int dPhi, int flags);
EXTERNC uint8 otsuThreshold(const unsigned char *image, const int size);
//EXTERNC int IBD_RefinePupil(SSegmentationResult *Result, const uint8* img, int H, int W, int flags);

EXTERNC int CircBypass(int* radPath, uint8* priceMap);