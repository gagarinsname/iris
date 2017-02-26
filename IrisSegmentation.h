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

//#include "ipl.h"
#include "unitypes.h"
#include "Geometry.h"

#define CONTOUR_QUANTIZATION 360

#define VOTING_LINES (0)
#define ANGDELT 17
#define RADDELT 2
#define MINRAD 20


typedef struct SEdgePnt
{
	int32 dir;
	int16 x;
	int16 y;
	int16 gx;
	int16 gy;
} SEdgePnt;

typedef struct SIrisData
{
	SCircleData *sCircle;
	double quality;
} SIrisData;

typedef struct SPupilData
{
	SCircleData *sCircle;
	SCircleData *sRCircle;
	sPoint* sContour;
	double quality;
} SPupilData;

typedef struct SSegmentationResult
{
	SPupilData *PupilData;
	SIrisData *IrisData;
	uint8* NormalizedIris;
} SSegmentationResult;


EXTERNC int IS_Initialize(SSegmentationResult *sResult);
EXTERNC int IrisSegmentation(SSegmentationResult* Result, //Segmentation result
	uint8 *imgInput, int H, int W, //image and its dimensions
	int angle, int flags);//stuff gonna be deleted
EXTERNC int IS_Deinitialize(SSegmentationResult *sResult);

EXTERNC int IS_ApproxBothBorders(SSegmentationResult *Result, uint8 *imgInput, int H, int W, int16 angle, int flags);

EXTERNC int pupilCircSP(int* radPath, uint8 *Map, int H, int W, int pRad);
EXTERNC int IS_RefinePupil(SSegmentationResult sResult, uint8 *img, int H, int W, int flags);