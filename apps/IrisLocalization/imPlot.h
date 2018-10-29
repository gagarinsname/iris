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
#include "unitypes.h"
#include "BorderDetection.h"
#include "stddefs.h"

// draw the line of the given color in grayscale image
EXTERNC void DRAW_Line(
	int* mat,	// image
	int H,		// image size
	int W,		// image size
	int x1,		// start point (x1, y1)
	int y1,		// start point (x1, y1)
	int x2,		// final point (x2, y2)
	int y2);	// final points (x2, y2)

// draw the circle of given color in grayscale image
EXTERNC void DRAW_CircleInGray(
	unsigned char* im,      // image
	int W,       // image size
	int H,       // image size
	SCircleData* sCircle, // circle structure
	uint8 color);    // color

void DRAW_2DSequenceInGray(
	uint8* im,		// image
	int W,			//image size
	int H,			//image size
	sPnt* sSeq,	//sequence points
	int N,			//Number of elements in the sequence
	uint8 color);	//color