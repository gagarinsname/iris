/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  25 Feb 2017                                            */
/*      Modified: many times                                             */
/*      Purpose:  2D data processing							         */
/*      Authors:                                                         */
/*        Yuriy Efimov                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

//#include "ipl.h"
#include "unitypes.h"
#include <math.h>
#pragma once

typedef struct sPoint
{
	int16 x, y;
} sPoint;

typedef struct SCircleData
{
	int32 xc;   // X coordinate
	int32 yc;   // Y coordinate
	int32 r;    // radius (elliptic big axis)
} SCircleData;
