/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  25 Feb 2017                                            */
/*      Modified: many times                                             */
/*      Purpose:  2D Geometry data processing							 */
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
	int16 xc;   // X coordinate
	int16 yc;   // Y coordinate
	int16 r;    // radius (elliptic big axis)
} SCircleData;
