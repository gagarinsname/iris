/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  26 December 2005                                       */
/*      Revision: 1.0.00                                                 */
/*      Purpose:  Image processing library                               */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#ifndef __OWNIPL_H__
#define __OWNIPL_H__

// error codes
#include "ipl.h"

/* image processing common error codes --------------------------------*/
/* queue overflow */
#define ERR_IPL_QUEUE_OVERFLOW ENHANCED_ERROR_CODE(0x202)
/* too many objects found - more specific than insufficient buffer */
#define ERR_IPL_TOO_MANY_OBJECTS ENHANCED_ERROR_CODE(0x203)
/* failed to restore gapped data - inadequate data on boundaries */
#define ERR_IPL_BAD_BOUNDARY_DATA ENHANCED_ERROR_CODE(0x206)
/* unexpressive data */
#define ERR_IPL_UNEXPRESSIVE_DATA ENHANCED_ERROR_CODE(0x207)
/* no objects found */
#define ERR_IPL_NO_OBJECTS_FOUND ENHANCED_ERROR_CODE(0x208)
/* sparse points */
#define ERR_IPL_SPARSE_POINTS ENHANCED_ERROR_CODE(0x209)

/* Calibration errors ------------------------------------------------*/
/* inadequate image for brightness calibration: too dark */
#define ERR_CLB_TOO_DARK_IMAGE ENHANCED_ERROR_CODE(0x241)
/* brightness calibration not needed: image is uniform */
#define ERR_CLB_BRTCALIB_NOT_NEEDED ENHANCED_ERROR_CODE(0x242)
/* inadequate image for brightness calibration: more than one big bright blobs */
#define ERR_CLB_MANY_BLOBS ENHANCED_ERROR_CODE(0x243)
/* black level mean square deviation is too big */
#define ERR_CLB_BLEVMSD_TOO_BIG ENHANCED_ERROR_CODE(0x244)
/* inadequate image for spatial calibration: not a square grid mesh */
#define ERR_CLB_NOT_SQUARE_GRID_MESH ENHANCED_ERROR_CODE(0x245)
/* inadequate image for spatial calibration: distortion center not in mesh */
#define ERR_CLB_CENTER_NOT_IN_MESH ENHANCED_ERROR_CODE(0x246)
/* inadequate image for spatial calibration: non convex cell found */
#define ERR_CLB_NON_CONVEX_CELL ENHANCED_ERROR_CODE(0x247)
/* inadequate image for spatial calibration: non-barrel distortion */
#define ERR_CLB_NON_BARREL ENHANCED_ERROR_CODE(0x248)
/* too big gaps in calibrating distances */
#define ERR_CLB_TOO_BIG_GAPS ENHANCED_ERROR_CODE(0x249)
/* low contrast */
#define ERR_CLB_LOW_CONTRAST ENHANCED_ERROR_CODE(0x24a)
/* too many edges in a place of cross */
#define ERR_CLB_TOO_MANY_EDGES ENHANCED_ERROR_CODE(0x24b)
/* too few edges in a place of cross */
#define ERR_CLB_TOO_FEW_EDGES ENHANCED_ERROR_CODE(0x24c)
/* not straight line cross */
#define ERR_CLB_NOT_STRAIGHT_LINE_CROSS ENHANCED_ERROR_CODE(0x24d)
/* too sparse mesh */
#define ERR_CLB_TOO_SPARSE_MESH ENHANCED_ERROR_CODE(0x24e)

// GRL (and not only) image and filter size limitations
#define PROCESSED_IMAGE_MAX_WIDTH 10240
#define PROCESSED_IMAGE_MAX_HEIGHT 10240
#define GRL_FILTER_MAX_HALFWIDTH 20
#define GRL_FILTER_MAX_HALFHEIGHT 20
#define GRL_FILTER_MAX_WIDTH (2*GRL_FILTER_MAX_HALFWIDTH+1)
#define GRL_FILTER_MAX_HEIGHT (2*GRL_FILTER_MAX_HALFHEIGHT+1)

void ownIPL_FillBorders1B(
  unsigned char* im,
  int xs,
  int ys);

void ownIPL_FillBorders1B_int(
  int* im,
  int xs,
  int ys);

void ownIPL_FillBordersXB(
  unsigned char* dst,
  int xs,
  int ys,
  int kehw,
  int kehh);

void ownIPL_FillBordersXB_int(
  int* dst,
  int xs,
  int ys,
  int kehw,
  int kehh);

#endif //__OWNIPL_H__
