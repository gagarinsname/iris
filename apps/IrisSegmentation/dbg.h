/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  4 April 2001                                           */
/*      Modified: 19 July 2010 - revised for CodCon                      */
/*      Revision: 2.0.00                                                 */
/*      Purpose:  debug facilitating stuff                               */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#ifndef __dbg_h__
#define __dbg_h__

#include <stdio.h>
#include "errorcodes.h"
#include "stddefs.h"
#include "bdl.h"

MIA_RESULT_CODE DBGL_GetVersionInfo(
  char* pcVersionString,  // OUT: place reserved for asciiz string
  int pnStringSize);      // IN:  string length

//== image pumps ==========================================================
// make 3Bpp from 1Bpp image
EXTERNC void DBGL_PUMP_Make3Bpp(
  unsigned char* dst,       // destination
  unsigned int dst_str,    // stride in bytes
  const unsigned char* src, // source
  unsigned int src_str,    // stride in bytes
  unsigned int w,          // ROI size
  unsigned int h);         //

// same as PumpMake3Bpp plus integral up-scaling
EXTERNC void DBGL_PUMP_Make3Bpp_Ext(
  unsigned char* dst,       // destination
  unsigned int dst_str,    // stride in bytes
  const unsigned char* src, // source
  unsigned int src_str,    // stride in bytes
  unsigned int w,          // ROI size
  unsigned int h,          //
  double scale);

// fill rgb24 image with given color
EXTERNC void DBGL_PUMP_FillRGB24(
  unsigned char* im,        // IN:  image
  int str,                  // IN:  stride in bytes
  int xs,                   // IN:  image size
  int ys,                   //
  int color);               // IN:  color

//== straight line figures ================================================
EXTERNC void DBGL_DRAW_LineInGray(
  unsigned char* im,    // image
  int str,    // stride in bytes
  int xs,     // image size in pixels
  int ys,     //
  unsigned int col,    // color
  int x0,     // line begin
  int y0,     //
  int x1,     // line end
  int y1);    //

EXTERNC void DBGL_DRAW_LineInRGB(
  unsigned char* im,    // image
  int str,    // stride in bytes
  int xs,     // image size in pixels
  int ys,     //
  unsigned int col,    // color as rgbquad
  int x0,     // line begin
  int y0,     //
  int x1,     // line end
  int y1);    //

EXTERNC void DBGL_DRAW_RectInGray(
  unsigned char* im,    // image
  int str,    // stride in bytes
  int xs,     // image size in pixels
  int ys,     //
  unsigned int col,    // color as rgbquad
  int x,      // rect corner
  int y,      //
  int w,      // rect size
  int h);     //

EXTERNC void DBGL_DRAW_RectInRGB(
  unsigned char* im,    // image
  int str,    // stride in bytes
  int xs,     // image size in pixels
  int ys,     //
  unsigned int col,    // color as rgbquad
  int x,      // rect corner
  int y,      //
  int w,      // rect size
  int h);     //

// draw cross sign of given color in grayscale (8 bpP) image
EXTERNC void DBGL_DRAW_CrossInGray(
  unsigned char* im,        // IN: image to be marked
  int xs,                   // IN: image size
  int ys,                   //
  int x,                    // IN: center of cross sign
  int y,                    //
  int sz,                   // IN: half size of the sign
  unsigned char color);     // color of the sign

EXTERNC void DBGL_DRAW_CrossInRGB24(
  unsigned char* im,        // IN: image to be marked
  int xs,                   // IN: image size
  int ys,                   //
  int x,                    // IN: center of cross sign
  int y,                    //
  int sz,                   // IN: half size of the sign
  unsigned int color);     // color of the sign

EXTERNC void DBGL_DRAW_BarInRGB(
  unsigned char* im,    // image
  int str,    // stride in bytes
  int xs,     // image size in pixels
  int ys,     //
  int xc,      // center
  int yc,      //
  int w,      // half sizes
  int h,      //
  unsigned int col);   // color as rgbquad

// draw plus sign of given color in grayscale (8 bpP) image
EXTERNC void DBGL_DRAW_PlusInGray(
  unsigned char* im,        // IN: image to be marked
  int xs,                   // IN: image size
  int ys,                   //
  int x,                    // IN: center of plus sign
  int y,                    //
  int sz,                   // IN: half size of the sign
  unsigned char color);     // color of the sign

// draw plus sign of given color in color (RGB24) image
EXTERNC void DBGL_DRAW_PlusInRGB24(
  unsigned char* im,        // IN: image to be marked
  int xs,                   // IN: image size
  int ys,                   //
  int x,                    // IN: center of plus sign
  int y,                    //
  int sz,                   // IN: half size of the sign
  int color);               // color of the sign

// draw plus sign of given color in color (RGB24) image
EXTERNC void DBGL_DRAW_PlusInRGBA(
  uint32* im,        // IN: image to be marked
  int xs,                   // IN: image size
  int ys,                   //
  int x,                    // IN: center of plus sign
  int y,                    //
  int sz,                   // IN: half size of the sign
  uint32 color);            // color of the sign

// draw plus sign of given color in color (RGB24) image
EXTERNC void DBGL_DRAW_PlusInScaledRGB24(
  unsigned char* im,        // IN: image to be marked
  int xs,                   // IN: image size
  int ys,                   //
  int x,                    // IN: center of plus sign
  int y,                    //
  int sz,                   // IN: half size of the sign
  int color,                // IN: color of the sign
  double dScale);              // IN: 'im' is supposed to be upscaled by this factor

//== round figures ===============================================================
EXTERNC void DBGL_DRAW_CircleInGray(
  unsigned char* im,      // image
  int xs,       // image size
  int ys,       // image size
  int xc,       // circle center
  int yc,       // circle center
  int r,        // circle radius
  unsigned char color);   // color

// draw the circle of xoir color in grayscale image
EXTERNC void DBGL_DRAW_CircleXORInGray(
  unsigned char* im,      // image
  int xs,       // image size
  int ys,       // image size
  int xc,       // circle center
  int yc,       // circle center
  int r);        // circle radius

EXTERNC void DBGL_DRAW_FillCircleInGray(
  unsigned char* im,      // image
  int xs,       // image size
  int ys,       // image size
  int xc,       // circle center
  int yc,       // circle center
  int r,        // circle radius
  unsigned char color);   // color

// draw the circle of given color on RGB24 image
EXTERNC void DBGL_DRAW_CircleInRGB(
  unsigned char* buf,     // image
  int str,      // stride in bytes
  int xs,       // image size
  int ys,       // image size
  int xc,       // circle center
  int yc,       // circle center
  int r,        // circle radius
  unsigned int color);   // RGBQUAD color

// draw the circle of given color on RGBA image
EXTERNC void DBGL_DRAW_CircleInRGBA(
  uint32* buf,     // image
  int str,      // stride in elements
  int xs,       // image size
  int ys,       // image size
  int xc,       // circle center
  int yc,       // circle center
  int r,        // circle radius
  uint32 color);   // RGBQUAD color

// draw the circle of given color on RGB24 image
EXTERNC void DBGL_DRAW_ThickCircleInRGB(
  unsigned char* buf,     // image
  int str,      // stride in bytes
  int xs,       // image size
  int ys,       // image size
  int xc,       // circle center
  int yc,       // circle center
  int r,        // circle radius
  unsigned int color,    // RGBQUAD color
  int thickness);

// draw ellipse on RGB24 image
EXTERNC void DBGL_DRAW_EllipseInRGB24(
  unsigned char* dst,   // destination image (RGB)
  int xs,     // image size
  int ys,     //
  double xc,      // center position
  double yc,      //
  double a,      // half-size
  double b,      //
  double f,     // tilt
  int color); // color in RGBQUAD

//== different figures ============================================================
// draw the contour given in polar coordinates as a function 'ro(phi)'
EXTERNC void DBGL_DRAW_PolarContourInRGB(
  unsigned char* im,    // destination image
  int xs,               // image size
  int ys,               //
  int xc,               // center of contour
  int yc,               //
  const int* rads,      // ro(phi) table
  int nPoints,          // number of points
  double dBegAngle,     // starting angle
  double dEndAngle,     // ending angle
  int color,            // color in image
  int isPerversed);     // =1 - count phi clockwise

// draw the contour given in polar coordinates as a function 'ro(phi)'
EXTERNC void DBGL_DRAW_PolarContourWithQualityInRGB(
  unsigned char* im,    // destination image
  int xs,               // image size
  int ys,               //
  int xc,               // center of contour
  int yc,               //
  const int* rads,      // ro(phi) table
  const int* quals,     // Q(phi) table in range [0-100]
  int nPoints,          // number of points
  double dBegAngle,     // starting angle
  double dEndAngle,     // ending angle
  int colorBest,        // color of best quality point in image
  int colorWorst,       // color of worst quality point in image
  int isPerversed);     // =1 - count phi clockwise

// draw the contour given by Cartesian points
EXTERNC void DBGL_DRAW_CartesianContourWithQualityInRGB(
  unsigned char* im,    // destination image
  int xs,               // image size
  int ys,               //
  const int* Xs,        // x table
  const int* Ys,        // y table
  const int* quals,     // Q(phi) table in range [0-100]
  int nPoints,          // number of points
  int colorBest,        // color of best quality point in image
  int colorWorst);      // color of worst quality point in image

// draw curve ax2+bxy+cy2+dx+ey+f=0
//            0   1   2   3  4  5
EXTERNC void DBGL_DRAW_GraphCurveOrder2InRGB(
  unsigned char* im,    // image
  int str,              // stride in bytes
  int xs,               // width
  int ys,               // height
  int color,            // color of the curve
  const double* coefs); // coefficients a,b,c,d,e,f

//== histograms =================================================================
EXTERNC void DBGL_DRAW_HistogramInGray(
  unsigned char* im,
  int xs,
  int ys,
  const int* pHist,
  int hLen,
  int bNorm);

EXTERNC void DBGL_DRAW_HistogramInRGBext(
  unsigned char* im,
  int xs,
  int ys,
  int str,
  const int* pHist,
  int hLen,
  int mode);

EXTERNC void DBGL_DRAW_HistogramInRGBext_uchar(
  unsigned char* im,
  int xs,
  int ys,
  int str,
  const unsigned char* pHist,
  int hLen,
  int mode);

EXTERNC void DBGL_DRAW_HistogramInGrayExt(
  unsigned char* im,  // IN:  buffer for image of the histogram
  int xs,             // IN:  histogram image dimensions
  int ys,             //
  int str,            // IN:  stride
  const int* pHist,   // IN:  histogram data
  int hLen,           // IN:  histogram length
  int mode);          // IN:  1 -normalize to max data value

EXTERNC void DBGL_DRAW_HistogramInGray_Double(
  unsigned char* im,    // IN:  buffer for image of the histogram
  int xs,               // IN:  histogram image dimensions
  int ys,               //
  int str,              // IN:  stride
  const double* pHist,  // IN:  histogram data
  int hLen);            // IN:  histogram length

EXTERNC void DBGL_DRAW_ColorHistogram(
  unsigned char* im,
  int xs,
  int ys,
  int str,
  const int* pHist,
  const int* pColors,
  int hLen,
  int mode);

EXTERNC void DBGL_DRAW_StrobesInRGB(
  unsigned char* im,
  int xs,
  int ys,
  int str,
  const int* pHist,
  int hLen,
  const int* pStrobeList,
  int nStrobeNum,
  int color);

//== text in image ================================================================
EXTERNC void DBGL_TYPE_CharInGray(
  unsigned char *im,  // IN:  image to draw in
  int xs,             // IN:  image size
  int ys,             // IN:  image size
  char ch,            // IN:  char to be printed
  int xb,             // IN:  position
  int yb,             // IN:  position
  int fgColor,        // IN:  foreground color (-1 for transparent)
  int bgColor);       // IN:  background color (-1 for transparent)

EXTERNC void DBGL_TYPE_CharInRGB(
  unsigned char *im,  // IN:  image to draw in
  int xs,             // IN:  image size
  int ys,             // IN:  image size
  char ch,            // IN:  char to be printed
  int xb,             // IN:  position
  int yb,             // IN:  position
  int fgColor,        // IN:  foreground color (-1 for transparent)
  int bgColor);       // IN:  background color (-1 for transparent)

EXTERNC void DBGL_TYPE_CharInRGBA(
  uint32 *im,         // IN:  image to draw in
  int xs,             // IN:  image size
  int ys,             // IN:  image size
  char ch,            // IN:  char to be printed
  int xb,             // IN:  position
  int yb,             // IN:  position
  int fgColor,        // IN:  foreground color (-1 for transparent)
  int bgColor);       // IN:  background color (-1 for transparent)

EXTERNC void DBGL_TYPE_TextInGray(
  unsigned char *im,  // IN:  image to draw in
  int xs,             // IN:  image size
  int ys,             // IN:  image size
  const char *str,    // IN:  string to be printed
  int xb,             // IN:  start position
  int yb,             // IN:  start position
  int fgColor,        // IN:  foreground color (-1 for transparent)
  int bgColor);       // IN:  background color (-1 for transparent)

EXTERNC void DBGL_TYPE_TextInRGB(
  unsigned char *im,  // IN:  image to draw in
  int xs,             // IN:  image size
  int ys,             // IN:  image size
  const char *str,    // IN:  string to be printed
  int xb,             // IN:  start position
  int yb,             // IN:  start position
  int fgColor,        // IN:  foreground color (-1 for transparent)
  int bgColor);       // IN:  background color (-1 for transparent)

EXTERNC void DBGL_TYPE_TextInRGBA(
  uint32* im,         // IN:  image to draw in
  int xs,             // IN:  image size
  int ys,             // IN:  image size
  const char *str,    // IN:  string to be printed
  int xb,             // IN:  start position
  int yb,             // IN:  start position
  int fgColor,        // IN:  foreground color (-1 for transparent)
  int bgColor);       // IN:  background color (-1 for transparent)

//== saving images to files =====================================================

EXTERNC MIA_RESULT_CODE DBGL_FILE_SaveInt8Image(
  const int8* im,
  int str,
  int xs,
  int ys,
  const char* filename);

EXTERNC MIA_RESULT_CODE DBGL_FILE_SaveUint8Image(
  const uint8* im,
  int str,
  int xs,
  int ys,
  const char* filename);

EXTERNC MIA_RESULT_CODE DBGL_FILE_SaveInt16Image(
  const int16* im,
  int str,
  int xs,
  int ys,
  const char* filename);

EXTERNC MIA_RESULT_CODE DBGL_FILE_SaveUint16Image(
  const uint16* im,
  int str,              // stride in elements
  int xs,
  int ys,
  const char* filename);

EXTERNC MIA_RESULT_CODE DBGL_FILE_SaveInt32Image(
  const int32* im,
  int str,
  int xs,
  int ys,
  const char* filename);

EXTERNC MIA_RESULT_CODE DBGL_FILE_SaveUint32Image(
  const uint32* im,
  int str,
  int xs,
  int ys,
  const char* filename);

EXTERNC MIA_RESULT_CODE DBGL_FILE_SaveDoubleImage(
  const double* im,
  int str,                // stride in elements
  int xs,
  int ys,
  const char* filename);

EXTERNC MIA_RESULT_CODE DBGL_FILE_SaveComplexImage(
  const double* im,
  int str,                // stride in elements
  int xs,
  int ys,
  const char* filename);

EXTERNC MIA_RESULT_CODE DBGL_FILE_SaveRGB8Image(
  const uint8* im,
  int str,              // IN:  stride in elements (3 bytes)
  int xs,
  int ys,
  const char* filename);

EXTERNC MIA_RESULT_CODE DBGL_FILE_SaveRGBA8Image(
  const uint32* im,
  int str,              // IN:  stride in elements (4 bytes)
  int xs,
  int ys,
  const char* filename);

//==================== parameter file I/O API ======================

typedef struct
{
  int nUpX1;
  int nUpY1;
  int nUpX2;
  int nUpY2;
  int nDnX1;
  int nDnY1;
  int nDnX2;
  int nDnY2;
  SPupilInfo Pi;
  SIrisInfo  Ii;
  int nPoints;
  short asPupilRads[DET_BORD_RADS];
//  int asPupilQuals[DET_BORD_RADS];  // qualities in scale 0-100
  short asIrisRads [DET_BORD_RADS];
} SOcclusionInfo;

typedef struct
{
  int nTotalQuality;          // total quality. negative - unacceptable
  int nGazeDirectionQuality;  // quality by flash decentration
  int nPupilSharpnessQuality; // quality by pupil border sharpness
  int nJLSInfoQuality;        // informativity by compression
  int nJLSOcclusionQuality;   // occlusion by compression
  int nShortCirquitQuality;   // 'resistance' by short-circuit algorithm
  int nHoleQuality;           // presence of bi radial holes
  int nPointEnergyQuality;    // energy of points stretching from circle
  int nPointMassQuality;      // share of points participating in interpol
} SQualityMeasures;

typedef struct
{
  char acFilename[FILENAME_MAX];
  int  nImageIndex;
  //SFlashInfo  sFI;
  SCenterInfo sCI;
  SBazRadInfo sBI;
  SPupilInfo  sPI;
  SIrisInfo   sII;
  SOcclusionInfo sOI;
  SQualityMeasures sQM;
  //int nSharpness;
  //int resolution;
//  short anPupilRadii[DET_BORD_RADS];
//  short anIrisRadii[DET_BORD_RADS];
} SParamsTable;

// save parameters to text file
// Format of: 
//  [name] [3 - center] [7 - bazrad] [4 - pupil] [4 - iris] [4 - occlusion angles]
EXTERNC MIA_RESULT_CODE PARAMS_SaveCommonTxt(
  SParamsTable* psPT,     // IN:  parameters table
  int nEntrs,             // IN:  number of entries
  const char* params);    // IN:  file name

// Load data from parameters file from various formats.
// Differentiating between formats is done by header.
EXTERNC MIA_RESULT_CODE PARAMS_LoadCommonTxt(
  SParamsTable** ppPT,    // OUT: parameters table
  int* pnEntrs,           // OUT: number of entries
  const char* params);    // IN:  file name

// save refined border data to text file
// Format: 
//   [name] [pupil x,y,n] [iris x,y,n] N*[pupil r] N*[iris r]
EXTERNC MIA_RESULT_CODE PARAMS_SaveDetailedBordersTxt(
  SParamsTable* psPT,     // IN/OUT: parameters table
  int nEntrs,             // IN: number of entries
  const char* params);    // IN:  file name

// load refined border data from text file
// Format: 
//   [name] [pupil x,y,n] [iris x,y,n] N*[pupil r] N*[iris r]
EXTERNC MIA_RESULT_CODE PARAMS_LoadDetailedBordersTxt(
  SParamsTable** ppPT,    // OUT: parameters table
  int* pnEntrs,           // OUT: number of entries
  const char* params);    // IN:  file name

EXTERNC MIA_RESULT_CODE PARAMS_InitRecord(
  SParamsTable* psPT);

typedef struct
{
  char name[FILENAME_MAX];
  int REx,REy,LEx,LEy,Nx,Ny,Mx,My;
} SFaceFileData;

EXTERNC MIA_RESULT_CODE FACE_LoadParams(
  SFaceFileData** ppsFD,
  int* pnItems,
  const char* filename);

#endif //__dbg_h__
