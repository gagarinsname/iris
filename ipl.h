/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  15 March 2001                                          */
/*      Modified: 15 February 2006 - grl removed                         */
/*      Modified: 19 July 2010 - revised for CodCon                      */
/*      Revision: 3.0.00                                                 */
/*      Purpose:  Image processing library                               */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#ifndef __ipl_h__
#define __ipl_h__

// needed for correct handling of floats
#include <math.h>

// common headers
#include "stddefs.h"
#include "errorcodes.h"

//== structures =============================================================
// point with attributes in image
typedef struct
{
  double x;                     // position of edge in image
  double y;                     //
  int isDefined:1;    // is the point detected
  int is13Quadrant:1; // quadrant of node
}	SDot;

typedef enum
{
  IPL_BIT_ERODE,
  IPL_BIT_DILATE,
  IPL_BIT_SUBTRACT,
  IPL_BIT_AND
} IPL_BIT_OPERATION;

// result output for FindMax function
typedef struct
{
  int value;
  int x;
  int y;
} SValuePosition;

// == version ====================================================================
EXTERNC MIA_RESULT_CODE IPL_GetVersionInfo(
  char* pcVersionString,  // OUT: place reserved for asciiz string
  int pnStringSize);      // IN:  string length

//== calibration and adjustment functions ========================================
// calculate mean and stddev of chessboard cell pixel size
EXTERNC MIA_RESULT_CODE IPL_ADJ_CalculatePixDist(
  double* pdist_x,    // OUT: pixel size along x
  double* pdist_y,    // OUT: pixel size along y
  double* pdev,       // OUT: deviation
  const SDot* dots,   // IN:  array of dots (unsorted, unaligned)
  int num,  // IN:  size of array
  int xs,   // IN:  image size
  int ys);  // IN:

// combine transformations given by XY-maps
// C[x] = (B * A)[x] = B[A[x]]   (A is done first)
EXTERNC MIA_RESULT_CODE IPL_ADJ_CombineTransfomations(
  double* C_x,          // OUT: destination X map
  double* C_y,          // OUT: destination Y map
  const double* A_x,    // IN:  first source X
  const double* A_y,    // IN:  first source Y
  const double* B_x,    // IN:  second source X
  const double* B_y,    // IN:  second source Y
  int xs,     // IN:  transformation map size
  int ys);    //

// get grosses positions in image of chessboard
EXTERNC MIA_RESULT_CODE IPL_ADJ_GetChessCrosses2(
  SDot* pCrosses,           // OUT: crosses found
  int* num,       // IN/OUT: reserved space/number of actually found crosses
  const unsigned char* im,  // IN: source image
  int xs,         // IN: image width in pixels
  int ys,         // IN: image height in pixels
  unsigned char** bufs);    // temporary buffers

// subpixel refinement of points
EXTERNC MIA_RESULT_CODE IPL_ADJ_RefineChessCrosses(
  SDot* pCrosses,           // OUT: crosses
  int num,        // IN: number of crosses
  const unsigned char* im,  // IN: source image
  int xs,         // IN: image width in pixels
  int ys,         // IN: image height in pixels
  int hws,        // IN: recommened horizontal half-size of window
  int vws);       // IN: recommened vertical half-size of window
                            //     if zero, functions guesses by itself

// determine the central white blob automatically
EXTERNC MIA_RESULT_CODE IPL_ADJ_FindCentralBlob(
  SDot* Marker,               // OUT: central dot position
  const unsigned char* im,    // IN:  source image
  int xs,           // IN:  image size
  int ys);          //

// calculate the value showing deviation of mesh from ideally parallelogramic
EXTERNC MIA_RESULT_CODE IPL_ADJ_EvaluateParallelogramicMeshQuality(
  double* quality,        // OUT: evaluated quality
  double* uncertainty,    // OUT: uncertainty (due to holes in mesh)
  const SDot* mesh,       // IN:  mesh to be evaluated
  int nx,       // IN:  dimensions of the mesh
  int ny);      // IN:

// associate crosses in left and right images
EXTERNC MIA_RESULT_CODE IPL_ADJ_CorrespondCrosses2(
  SDot* pLeftCrosses,     // IN/OUT: array of crosses on left image/sorted array
  SDot* pRightCrosses,    // IN/OUT: array of crosses on right image/sorted array
                          //    first, associated items go with black quadrant 13
                          //    second, associated items go with black quadrant 24
                          //    items not associated are put in tail
  int* pnumlef, // IN/OUT: number of left crosses/number of black quad 13
  int* pnumrig, // IN/OUT: number of right crosses/number of black quad 24
  int xs,       // IN: image width in pixels
  int ys,       // IN: image height in pixels
  int xl,       // IN: correspondent point coordinates in images
  int yl,       // IN:    (these are used as a reference center)
  int xr,       // IN: 
  int yr);      // IN: 

//== bit image processing ======================================================
// bit image morphology
// WARNING: checking boundaries is on behalf of a caller !
EXTERNC MIA_RESULT_CODE IPL_BIT_Morphology(
  unsigned char* dst,             // destination
  unsigned int dststr,           // stride in bytes
  const unsigned char* src,       // source
  unsigned int srcstr,           // stride in bytes
  unsigned int xs,               // image size
  unsigned int ys,               //
  IPL_BIT_OPERATION operation,    // morphological operation
  unsigned int size);            // size of operation

EXTERNC MIA_RESULT_CODE IPL_BIT_ByteImageToBitImage(
  unsigned int* dst,        // dest image
  int dst_str,              // dest stride in DWORDs
  const unsigned char* src, // source image
  int src_str,              // source stride in bytes
  int width,                // width of source region to be converted
  int height);              // height

//== connected objects processing ==============================================

// primary integral (connected) object parameters
typedef struct
{
  int32 c;    // object number (color in painting functions)
  int32 xb;   // object starting point
  int32 yb;   //
  int32 M;    // mass Sum(1)
  int32 MX;   // first-order moments  : Sum(x)
  int32 MY;   //                      : Sum(y)
  int64 MXX;  // second order moments : Sum(x*x)
  int64 MXY;  //                      : Sum(x*y)
  int64 MYY;  //                      : Sum(y*y)
} SConObjPrimaryParams;

// data of straight chunk
typedef struct
{
  int n;
  int x;
  int y;
  int gx;
  int gy;
  int c;  // color (for debug)
} SChunkData;

typedef struct
{
  int S1;     // mass
  int id;     // ID (color)
  double xc;  // center
  double yc;  // 
  double phi; // angle
  double len; // length (diameter)
  double err; // error (MSD per point)
} CHUNK_Linear;

// get (a) linear (straight) chunks by (b) forest fire method by (c) coordinates
// if poor point is met, stop object enumeration by forest fire instantly
EXTERNC MIA_RESULT_CODE IPL_CHUNK_Linear_ForestFire_Coords_v1(
  CHUNK_Linear** ppsCL,     // OUT: array of segments
  int* pnCL,                // OUT: number of segments
  int32** ppnIdxs,          // OUT: list of point indexes (len = total number of enumerated points)
  int32** ppnCnts,          // OUT: list of index counts (len = number of enumerated objects+1)
  uint32** ppCoIm,          // OUT: 'colored' image of segments (32bpP)
  const unsigned char* im,  // IN:  source image
  const short* gxim,        // IN:  brightness gradients of source image
  const short* gyim,        //
  int xs,                   // IN:  source image size
  int ys,                   //
  int nMassThr);            // IN:  object mass threshold

// get (a) linear (straight) chunks by (b) forest fire method by (c) coordinates
// if poor point is met, only this point is omitted
EXTERNC MIA_RESULT_CODE IPL_CHUNK_Linear_ForestFire_Coords_v2(
  CHUNK_Linear** ppsCL,     // OUT: array of segments
  int* pnCL,                // OUT: number of segments
  int32** ppnIdxs,          // OUT: list of point indexes (len = total number of enumerated points)
  int32** ppnCnts,          // OUT: list of index counts (len = number of enumerated objects+1)
  uint32** ppCoIm,          // OUT: 'colored' image of segments (32bpP)
  const unsigned char* im,  // IN:  source image
  const short* gxim,        // IN:  brightness gradients of source image
  const short* gyim,        //
  int xs,                   // IN:  source image size
  int ys,                   //
  int nMassThr);            // IN:  object mass threshold

// get (a) linear (straight) chunks by (b) forest fire method by (c) coordinates
// if poor point is met, it is postponed to end of queue
EXTERNC MIA_RESULT_CODE IPL_CHUNK_Linear_ForestFire_Coords_v3(
  CHUNK_Linear** ppsCL,     // OUT: array of segments
  int* pnCL,                // OUT: number of segments
  int32** ppnIdxs,          // OUT: list of point indexes (len = total number of enumerated points)
  int32** ppnCnts,          // OUT: list of index counts (len = number of enumerated objects+1)
  uint32** ppCoIm,          // OUT: 'colored' image of segments (32bpP)
  const unsigned char* im,  // IN:  source image
  const short* gxim,        // IN:  brightness gradients of source image
  const short* gyim,        //
  int xs,                   // IN:  source image size
  int ys,                   //
  int nMassThr);            // IN:  object mass threshold

// get (a) linear (straight) chunks by (b) hierarchical clustering from (c) coordinates
EXTERNC MIA_RESULT_CODE IPL_CHUNK_Linear_HierClust_Coords(
  CHUNK_Linear** ppsCL,     // OUT: array of segments
  int* pnCL,                // OUT: number of segments
  int32** ppnIdxs,          // OUT: list of point indexes (len = total number of enumerated points)
  int32** ppnCnts,          // OUT: list of index counts (len = number of enumerated objects+1)
  uint32** ppCoIm,          // OUT: 'colored' image of segments (32bpP)
  const unsigned char* im,  // IN:  source image
  const short* gxim,        // IN:  brightness gradients of source image
  const short* gyim,        //
  int xs,                   // IN:  source image size
  int ys,                   //
  int nMassThr);            // IN:  object mass threshold

// data of circular chunk (arc-chunk ~ archunk)
typedef struct
{
  int n;
  int x;
  int y;
  int r;
  int c;  // color (for debug)
} SArchunkData;

EXTERNC MIA_RESULT_CODE IPL_CONOB_ChunkDetect_v1(
  SChunkData** ppsCDO,      // OUT: array of segments
  int* pnCDO,               // OUT: number of segments
  const unsigned char* im,  // IN:  source image
  const short* gxim,        // IN:  brightness gradients of source image
  const short* gyim,        //
  int xs,                   // IN:  source image size
  int ys,                   //
  const char* nambeg);      // IN:  for debug, set to NULL

EXTERNC MIA_RESULT_CODE IPL_CONOB_ArchunkDetect_v1(
  SArchunkData** ppsCDI,    // OUT: array of arcs
  int* pnCDI,               // OUT: number of segments
  const unsigned char* im,  // IN:  source image
  const short* gxim,        // IN:  brightness gradients of source image
  const short* gyim,        //
  int xs,                   // IN:  source image size
  int ys,                   //
  const char* nambeg);      // IN:  for debug, set to NULL

// Non-optimized enumeration function (operating with pixels)
MIA_RESULT_CODE IPL_CONOB_EnumerateConnectedObjects_pixeled(
  SConObjPrimaryParams* pCOPP,  // OUT: buffer for object parameters
  int* pObjnum,        // IN/OUT: number of objects allocated/found
  unsigned char* im,            // input image (completely zeroed inside the function)
  int str,             // image stride in elements
  int xs,              // image size
  int ys,              //
  int MassThr_low,     // lower object mass threshold (-1 - disable)
  int MassThr_high,    // higher object mass threshold (-1 - disable)
  int connect8);       // if(TRUE) connect 8 neibours else - 4

// enumeration function optimized by segments
EXTERNC MIA_RESULT_CODE IPL_CONOB_EnumerateConnectedObjects(
  SConObjPrimaryParams* pCOPP,  // OUT: buffer for object parameters
  int* pObjnum,        // IN/OUT: number of objects allocated/found
  unsigned char* im,   // input image (completely zeroed inside the function)
  int str,             // image stride in elements
  int xs,              // image size
  int ys,              //
  int MassThr_low,     // lower object mass threshold (-1 - disable)
  int MassThr_high,    // higher object mass threshold (-1 - disable)
  int connect8);       // if(TRUE) connect 8 neibours else - 4

// Non-optimized function
EXTERNC MIA_RESULT_CODE IPL_CONOB_FloodFill(
  unsigned char* im,              // input image (completely zeroed inside the function)
  unsigned int str,               // image stride in elements
  unsigned int xs,                // image size
  unsigned int ys,                //
  unsigned int xc,                // center to begin fill
  unsigned int yc,                //
  unsigned int connect8,          // if(TRUE) connect 8 neibours else - 4
  unsigned int color);            // filling color

// Non-optimized painting function
EXTERNC MIA_RESULT_CODE IPL_CONOB_PaintConnectedObjects(
  unsigned int* painting,        // OUT: painted buffer
  unsigned int dststr,           // stride in elements
  SConObjPrimaryParams* pCOPP, // OUT: buffer for object parameters
  unsigned int* pObjnum,         // IN:  size of this buffer
                                  // OUT: number of objects found
  unsigned char* im,              // input image (completely zeroed inside the function)
  unsigned int str,              // image stride in elements
  unsigned int xs,               // image size
  unsigned int ys,               //
  unsigned int MassThr,          // object mass threshold
  unsigned int connect8);        // if(TRUE) connect 8 neibours else - 4

// Non-optimized painting function
MIA_RESULT_CODE IPL_CONOB_PaintConnectedObjectsThresholded(
  unsigned int* painting,       // OUT: painted buffer
  unsigned int dststr,          // stride in elements
  SConObjPrimaryParams* pCOPP,  // OUT: buffer for object parameters
  unsigned int* pObjnum,        // IN/OUT:  size of buffer/number of objects found
  unsigned char* im,            // input image (completely zeroed inside the function)
  unsigned int str,             // image stride in elements
  unsigned int xs,              // image size
  unsigned int ys,              //
  unsigned int MassThr,         // object mass threshold
  unsigned int connect8,        // if(TRUE) connect 8 neibors else - 4
  uint8 ThrActivate,
  uint8 ThrContinue);

//== some drawing used for image rocessing ======================================
EXTERNC void IPL_DRAW_Line(
  unsigned char* im,    // image
  unsigned int str,    // stride in bytes
  unsigned int xs,     // image size in pixels
  unsigned int ys,     //
  unsigned int col,    // color
  unsigned int x0,     // line begin
  unsigned int y0,     //
  unsigned int x1,     // line end
  unsigned int y1);    //

EXTERNC MIA_RESULT_CODE IPL_DRAW_MaskFromContour(
  unsigned char** pMask,
  int* pX,
  int* pY,
  int* pW,
  int* pH,
  int xc,
  int yc,
  int Count,
  const int* Radii,
  int isClockWise);

EXTERNC MIA_RESULT_CODE IPL_DRAW_MaskFromTwoContours(
  unsigned char* mask,  // size = xs*ys
  int xs,
  int ys,
  int isClockwise,
  int isPupilIncl,
  //const SEyeRefinedData* pERD,
  int xc_include,
  int yc_include,
  int Count_include,
  const int* Radii_include,
  int xc_exclude,
  int yc_exclude,
  int Count_exclude,
  const int* Radii_exclude);

//== ellipse detection ======================================
typedef struct
{
  double ellpar[5];
  double Qfill;     // fill quality
  double Qdir;      // direction quality
  double Qelong;    // elongation quality
  int binval;       // binarization value
  int multi;        // multiplicity factor
} SEllipseExt;

typedef struct
{
  double QfillMin;      // minimum acceptable fill
  double QelongMin;     // minimum acceptable elongation
  int nMassMin;         // minimum acceptable object mass
  int nMassMax;         // maximum acceptable object mass
  int aThresholds[256]; // list of threshold values
  int nThresholds;      // number of thresholds
} SEllDetPar;

// get the best ellipse
EXTERNC MIA_RESULT_CODE IPL_ELL_Detect_v1(
  double* pdElls,   // OUT: 5 ellipse parameters - x,y,a,b,f
  const uint8* im,  // IN:  source image
  int xs,           // IN:  image size
  int ys);

// detect multiple ellipses
EXTERNC MIA_RESULT_CODE IPL_ELL_DetectMulti_v1(
  SEllipseExt** ppsElls,  // OUT: set of detected ellipses with qualities
  int* pnElls,            // OUT: number of detected ellipses
  const uint8* im,        // IN:  source image
  int xs,                 // IN:  image size
  int ys);

// detect multiple ellipses
EXTERNC MIA_RESULT_CODE IPL_ELL_DetectMultiExt_v1(
  SEllipseExt** ppsElls,  // OUT: set of detected ellipses with qualities
  int* pnElls,            // OUT: number of detected ellipses
  const uint8* im,        // IN:  source image
  int xs,                 // IN:  image size
  int ys,
  const SEllDetPar* pEp); // IN:  parameters of ellipse selection

/*---------------------------------------------------------------------------*/
/* Purpose:  This routine replaces DATA by its one-dimensional discrete      */
/*           transform if ISIGN=1 or replaces DATA by its inverse transform  */
/*           if ISIGN=-1.  DATA is a complex array of length NN which is     */
/*           input as a real array of length 2*NN.  NN must be an integer    */
/*           power of 2 or the routine will abort and return INVALID.        */
/*                                                                           */
/* Note:     Because this code was adapted from a FORTRAN library, the       */
/*           data array is 1-indexed.  In other words, the first element     */
/*           of the array is assumed to be in data[1].  Because C is zero    */
/*           indexed, the first element of the array is in data[0].  Hence,  */
/*           we must subtract 1 from the data address at the start of this   */
/*           routine so references to data[1] will really access data[0].    */
/*---------------------------------------------------------------------------*/
EXTERNC int IPL_FFT_1D(double *data, int nn, int isign);

/*---------------------------------------------------------------------------*/
/* Purpose:  This routine replaces DATA by its two-dimensional discrete      */
/*           transform if ISIGN=1 or replaces DATA by its inverse transform  */
/*           if ISIGN=-1.  DATA is a complex array with NN columns and MM    */
/*           rows.  NN and MM must both be integer powers of 2 or the        */
/*           routine will abort and return INVALID.                          */
/*---------------------------------------------------------------------------*/
EXTERNC int IPL_FFT_2D(double data[], int nn, int mm, int isign);

/*---------------------------------------------------------------------------*/
/* Purpose:  This routine replaces DATA by its three-dimensional discrete    */
/*           transform if ISIGN=1 or replaces DATA by its inverse transform  */
/*           if ISIGN=-1.  DATA is a complex array with NN columns and MM    */
/*           rows and OO slices.  NN, MM and OO must be integer powers of    */
/*           2 or the routine will abort and return INVALID.                 */
/*---------------------------------------------------------------------------*/
EXTERNC int IPL_FFT_3D(double data[], int nn, int mm, int oo, int isign);

/*---------------------------------------------------------------------------*/
/* Purpose:  This routine replaces DATA by its DIMC-dimensional discrete     */
/*           Fourier transform, if SIGN is input as 1.  DIMV[0..DIMC-1]      */
/*           is an integer array containing the lengths of each dimension    */
/*           (number of complex values) , which MUST be a power of 2.        */
/*           DATA is a real array of length twice the product of these       */
/*           lengths, in which the data is stored in a multidimensional      */
/*           complex array: real and imaginary parts of this array are in    */
/*           consecutive locations, and the rightmost index of the array     */
/*           increases most rapidly as one proceeds along DATA.  For a two-  */
/*           dimensional array, this is equivalent to storing the array      */
/*           by rows (the C norm).  If SIGN is -1, DATA is replaced by its   */
/*           inverse transform times the product of the lengths of all       */
/*           dimensions (in other words, divide by the number of pixels      */
/*           in the image to get the inverse transform).                     */
/*                                                                           */
/* Note:     Because this code was adapted from a FORTRAN library, the       */
/*           data array is 1-indexed.  In other words, the first element     */
/*           of the array is assumed to be in data[1].  Because C is zero    */
/*           indexed, the first element of the array is in data[0].  Hence,  */
/*           the address of data[-1] must be passed to this routine to       */
/*           ensure that references to data[1] will really access data[0].   */
/*---------------------------------------------------------------------------*/
EXTERNC int IPL_FFT_ND(double data[], int nn[], int ndim, int isign);

//== filters ===================================================================

// anisotropic diffusion to enhance edges
EXTERNC MIA_RESULT_CODE IPL_FILT_AnisotropicDiffusion(
  double* dst,
  const double* src,
  int xs,
  int ys,
  double K,
  double lambda,
  int iterations,
  const char* nambeg);

// relaxation filter
EXTERNC MIA_RESULT_CODE IPL_FILT_RFilter(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  int kehw,
  int kehh);

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

EXTERNC void IPL_FILT_UniformBlur3x3_Int(
  int* dst,
  const int* src,
  int xs,
  int ys,
  int mode);  // 0 - normalized, 1- no normalization

// 5x5 Haussian blur
EXTERNC void IPL_FILT_HaussBlur5x5(
  unsigned char* dst,         // OUT: dst image
  const unsigned char* src,   // IN:  src image
  int xs,                     // IN:  image size
  int ys);                    //

EXTERNC void IPL_FILT_HaussBlur5x5_int32(
  int* dst,
  const int* src,
  int xs,
  int ys,
  int bNormalise);

// one dimentional haussian blurring 1x5
EXTERNC void IPL_FILT_HaussBlur1x5(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys);

// one dimentional haussian blurring 5x1
EXTERNC void IPL_FILT_HaussBlur5x1(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys);

// one dimentional vertical haussian blurring
EXTERNC void IPL_FILT_HaussBlur1x5(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys);

// one dimentional vertical haussian blurring
EXTERNC void IPL_FILT_ShiftedDiff1x3(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  int shift);                   // shift for difference: simulate gabor filter

// Blur an image with a Haussian filter.
EXTERNC MIA_RESULT_CODE IPL_FILT_GaussianSmooth_uint8(
  short *smoothedim,
  double *tempim,
  const unsigned char *image,
  int xs,
  int ys,
  double sigma);

// Blur an image with a gaussian filter
EXTERNC MIA_RESULT_CODE IPL_FILT_GaussianSmooth_ext_uint8(
  short *smoothedim,
  double *tempim,
  const unsigned char *image,
  int xs,
  int ys,
  double sigma_x,
  double sigma_y,
  double _BOOSTBLURFACTOR);

// Blur an image with a Haussian filter.
EXTERNC MIA_RESULT_CODE IPL_FILT_GaussianSmooth_double(
  double *smoothedim,
  double *tempim,
  const double *image,
  int xs,
  int ys,
  double sigma);

// calculation of mean in image
EXTERNC MIA_RESULT_CODE IPL_FILT_CalculateMean(
  unsigned char* dst,         // dst
  int dststr,                // dest stride
  const unsigned char* src,   // src
  int srcstr,   // source stride
  int xs,       // xs
  int ys,       // ys
  int x,        // ROI
  int y,        //
  int w,        //
  int h,        //
  int kew,      // kernel half size
  int keh);     //

// calculation of mean and dispersion on image
MIA_RESULT_CODE IPL_FILT_CalculateMeanDisp(
  unsigned char* mean,       // dst
  unsigned char* disp,       // dst
  int dststr,              // dest stride
  const unsigned char* src, // src
  int srcstr,   // source stride
  int xs,       // xs
  int ys,       // ys
  int x,        // ROI
  int y,        //
  int w,        //
  int h,        //
  int kew,      // kernel half size
  int keh,      //
  int dispmul); // dispersion multyplier

// calculate local variation
EXTERNC MIA_RESULT_CODE IPL_FILT_LocalVariation(
  unsigned char* dst,       // OUT: dst image
  const unsigned char* src, // IN:  source image
  unsigned int xs,          // IN:  image size
  unsigned int ys,          //
  unsigned int wx,          // IN:  window size
  unsigned int wy);         //

// calculate (val)-(min in local nbhood)
EXTERNC MIA_RESULT_CODE IPL_FILT_LocalMaxDrop(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  int wx,
  int wy);

// normalise globally (i.e. to fixed total brightness)
EXTERNC void IPL_FILT_GlobalNormalise(
  unsigned char* dst,         // destination
  int dststr,                 // stride
  const unsigned char* srcs,  // source (NULL for inplace)
  int srcstr,                 // stride
  int w,                      // ROI size
  int h,                      //
  unsigned char v);           // norm

// normalize locally (i.e. to fixed local brightness)
EXTERNC MIA_RESULT_CODE IPL_FILT_LocalNormalise(
  unsigned char* dst,         // dst
  int dststr,                // dest stride
  const unsigned char* src,   // src
  int srcstr,   // source stride
  int xs,       // xs
  int ys,       // ys
  int x,        // ROI
  int y,        //
  int w,        //
  int h,        //
  int kew,      // kernel half size
  int keh);     //

// detect chessboard nodes in image
EXTERNC MIA_RESULT_CODE IPL_FILT_ChessCrossesWithBlack13Quadrant(
  unsigned char* dst,         // destination image, stride=xs
  const unsigned char* er,    // source image (eroded)
  const unsigned char* di,    // source image (dilated)
  unsigned int srcstr,       // its stride
  unsigned int xs,           // image width
  unsigned int ys,           // image height
  unsigned int kew);         // kernel half-size

// detect chessboard nodes in image
EXTERNC MIA_RESULT_CODE IPL_FILT_ChessCrossesWithBlack24Quadrant(
  unsigned char* dst,         // destination image, stride=xs
  const unsigned char* er,    // source image (eroded)
  const unsigned char* di,    // source image (dilated)
  unsigned int srcstr,       // its stride
  unsigned int xs,           // image width
  unsigned int ys,           // image height
  unsigned int kew);         // kernel half-size

// 3x3 erosion
EXTERNC void IPL_FILT_Erode3x3(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys);

// 3x3 dilatation
EXTERNC void IPL_FILT_Dilate3x3(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys);

// '+'-shaped mask
EXTERNC void IPL_FILT_Erode3x3Cross(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys);

// '+'-shaped mask
EXTERNC void IPL_FILT_Dilate3x3Cross(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys);

EXTERNC void IPL_FILT_Erode5x1( //  along x
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys);

EXTERNC void IPL_FILT_Dilate5x1( //  along x
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys);

EXTERNC void IPL_FILT_Erode1x5( //  along y
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys);

EXTERNC void IPL_FILT_Dilate1x5( //  along y
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys);

// Dilated image should be binary, i.e. contain only two values: 
// zero and non-zero (and ALL non-zero values should be equal)
// 255 is used in dst as non-zero value. 
EXTERNC MIA_RESULT_CODE IPL_FILT_Dilate_Binary(
  unsigned char* dst,         // OUT: resulting image
  const unsigned char* src,   // IN:  original image
  int xs,                     // IN:  image size
  int ys,                     //
  int kehw,                   // IN:  kernel half-size
  int kehh,                   //
  void* pvBuf,                // IN:  temporary buffer
  int* pnLen);                // IN/OUT: bytes allocated/used

// 1D erosion for int type
EXTERNC void IPL_FILT_Erode1D_Cycled_Int(
  int* dst,         // OUT: destination array
  const int* src,   // IN:  source array
  int len,          // IN:  size of array
  int kehw);        // IN:  half aperture

// 1D dilation for int type
EXTERNC void IPL_FILT_Dilate1D_Cycled_Int(
  int* dst,         // OUT: destination array
  const int* src,   // IN:  source array
  int len,          // IN:  size of array
  int kehw);        // IN:  half aperture

// general median filter
EXTERNC void IPL_FILT_Median(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  int kehw,
  int kehh);

EXTERNC void IPL_FILT_Dilate(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  int kehw,
  int kehh);

EXTERNC void IPL_FILT_Erode(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  int kehw,
  int kehh);

// median in 3x3 window (by Knuth)
EXTERNC MIA_RESULT_CODE IPL_FILT_Median3x3(
  unsigned char* dst,           // destination
  const unsigned char* src,     // source
  int xs,               // size of image in elements
  int ys);

// median in 3x1 window for ints
EXTERNC void IPL_FILT_Median3x1_Int(
  int* dst,
  const int* src,
  int xs,
  int ys);

EXTERNC void IPL_FILT_mdn2(
  const unsigned char *bufin,
  unsigned char *bufout,
  int lstr,
  int apt_size);

EXTERNC float IPL_FILT_Mean_Disp(
  short *bufin,
  short *Midl,
  int lstr);

EXTERNC void IPL_FILT_mdnS256(
  short *bufin,
  short *bufout,
  int lstr,
  int apt_size);

// straight-thru median for sequence of int's
EXTERNC MIA_RESULT_CODE IPL_FILT_LineMedian_Int(
  int *dst,         // OUT: output array
  const int* src,   // IN:  input array
  int len,          // IN:  length of input
  int hw);          // IN:  half-width of median window

// straight-thru median for sequence of double's
EXTERNC MIA_RESULT_CODE IPL_FILT_LineMedian_Double(
  double *dst,        // OUT: output array
  const double* src,  // IN:  input array
  int len,            // IN:  length of input
  int hw);            // IN:  half-width of median window

// canny edge detection
EXTERNC MIA_RESULT_CODE IPL_FILT_Canny(
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

// edge drawing
EXTERNC MIA_RESULT_CODE IPL_FILT_EdgeDraw_v1(
  unsigned char *edge,        // OUT: canny mask
  const unsigned char *image, // IN:  source image
  short ** ppgx,              // OUT: map of horizontal gradients
  short ** ppgy,              // OUT: map of vertical gradients
  int xs,                     // IN:  image size
  int ys,                     // 
  void* pvBuf,                // IN:  temporary buffer
  int* pnLen,                 // IN/OUT: bytes allocated/used
  const char *nambeg);        // IN:  set this to NULL

EXTERNC void IPL_FILT_ShiftedDiff1x5(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  int shift);                 // shift for difference: simulate gabor filter

EXTERNC void IPL_FILT_ShiftedDiff5x1(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  int shift);                 // shift for difference: simulate gabor filter

//== geometry functions (changing images not-locally) ==========================
EXTERNC MIA_RESULT_CODE IPL_GEOM_GetPolarImage0(
  const unsigned char *BufImg,// IN:  eye image
  int     SizeX,        // IN: image size
  int     SizeY,        //
  int     SizeStride,        //
  unsigned char  *bufpol,   // IN-OUT: polar data block 
  int            SizePolarA,   // IN: number of lines in Polar template
  int            SizePolarR,   // IN: number of columns in Polar template
  int            SizePolarStride, // IN: number of columns in Polar template
  int            x,            //center coordinates
  int            y,
  int     BegRad,       // OUT: beginning radius polar data block 
  float   BegAng,       // IN:  rotation angle for frame
  void*   MemoryBuff,   // IN:  tmp memory buffer
  int*     MemoryBuffSize);

// using non-concentric pupil and iris circles
EXTERNC MIA_RESULT_CODE IPL_GEOM_GetPolarImage1(  // MIA 11.08.04
  unsigned char* dst, // OUT: polar data block 
  int phi,            // IN: number of lines in polar image
  int rho,            // IN: number of columns in polar image
  const unsigned char* src,// IN:  eye image
  int xs,             // IN: image size
  int ys,             //
  int xp,             // center coordinates of pupil (inner boundary)
  int yp,             //
  int rp,             // radius of pupil
  int xi,             // center coordinates of iris (outer boundary)
  int yi,             //
  int ri);            // radius of iris

// rotate and resize ROI with bi-linear interpolation (preserving orthogonality)
EXTERNC MIA_RESULT_CODE IPL_GEOM_RotateResizeBLI(
  unsigned char* dst,       // destination image
  int dststr,               // stride in bytes
  const unsigned char* src, // source image
  int srcstr,               // stride in bytes
  int xs,                   // source image size
  int ys,                   //
  float x,                  // source ROI
  float y,                  //
  float w,                  //
  float h,                  //
  float a,                  // angle of source ROI tilt
  int W,                    // destination ROI size
  int H);                   //

// rotate and resize ROI with bi-linear interpolation 
// using affine transform matrix (no orthogonality)
EXTERNC MIA_RESULT_CODE IPL_GEOM_AffineTransformBLI(
  unsigned char* dst,       // destination image
  unsigned int dststr,      // stride in bytes
  const unsigned char* src, // source image
  unsigned int srcstr,      // stride in bytes
  unsigned int xs,          // source image size
  unsigned int ys,          //
  double x,                 // corner of source ROI
  double y,                 //
  const double* A,          // affine transformation matrix
  unsigned int W,           // destination ROI size
  unsigned int H);          //

// rotate and resize ROI with bi-linear interpolation
EXTERNC MIA_RESULT_CODE IPL_GEOM_AffineTransformBLI_RGB(
  unsigned char* dst,         // destination image
  unsigned int dststr,        // stride in bytes
  const unsigned char* src,   // source image
  unsigned int srcstr,        // stride in bytes
  unsigned int xs,            // source image size
  unsigned int ys,            //
  double x,                   // source ROI
  double y,                   //
  const double* A,            // affine transformation matrix
  unsigned int W,             // destination ROI size
  unsigned int H);            //

// rotate and resize ROI with bi-linear interpolation
MIA_RESULT_CODE IPL_GEOM_AffineTransformBLI_RGBtoGRAY(
  unsigned char* dst,         // destination image
  unsigned int dststr,        // stride in bytes
  const unsigned char* src,   // source image
  unsigned int srcstr,        // stride in bytes
  unsigned int xs,            // source image size
  unsigned int ys,            //
  double x,                   // source ROI
  double y,                   //
  const double* A,            // affine transformation matrix
  unsigned int W,             // destination ROI size
  unsigned int H);            //

// YUV422 conversion
// rotate and resize ROI with bi-linear interpolation
MIA_RESULT_CODE IPL_GEOM_AffineTransformBLI_YUV422toGRAY(
	unsigned char* dst,         // destination image
	unsigned int dststr,        // stride in bytes
	const unsigned char* src,   // source image
	unsigned int srcstr,        // stride in bytes
	unsigned int xs,            // source image size
	unsigned int ys,            //
	double x,                   // source ROI
	double y,                   //
	const double* A,            // affine transformation matrix
	unsigned int W,             // destination ROI size
	unsigned int H);             //

// BAYER conversion
// rotate and resize ROI with bi-linear interpolation from Bayer RGGB
MIA_RESULT_CODE IPL_GEOM_AffineTransformBLI_BAYERtoGRAY(
	unsigned char* dst,         // destination image
	unsigned int dststr,        // stride in bytes
	const unsigned char* src,   // source image
	unsigned int srcstr,        // stride in bytes
	unsigned int xs,            // source image size
	unsigned int ys,            //
	double x,                   // source ROI
	double y,                   //
	const double* A,            // affine transformation matrix
	unsigned int W,             // destination ROI size
	unsigned int H);             //

// BAYER conversion
// R G R G R G
// G B G B G B
// rotate and resize ROI with nearest neighbour interpolation from Bayer RGGB
MIA_RESULT_CODE IPL_GEOM_AffineTransform_BAYERtoGRAYcrop(
	unsigned char* dst,         // destination image
	unsigned int dststr,        // stride in bytes
	const unsigned char* src,   // source image
	unsigned int srcstr,        // stride in bytes
	unsigned int xs,            // source image size
	unsigned int ys,            //
	double x,                   // source ROI
	double y,                   //
	const double* A,            // affine transformation matrix
	unsigned int W,             // destination ROI size
	unsigned int H);            //

// rotate and resize ROI with nearest-neighbor interpolation
EXTERNC MIA_RESULT_CODE IPL_GEOM_AffineTransformNNI_RGB(
  unsigned char* dst,         // destination image
  unsigned int dststr,        // stride in bytes
  const unsigned char* src,   // source image
  unsigned int srcstr,        // stride in bytes
  unsigned int xs,            // source image size
  unsigned int ys,            //
  double x,                   // source ROI
  double y,                   //
  const double* A,            // affine transformation matrix
  unsigned int W,             // destination ROI size
  unsigned int H);            //

// mirror image relative to vertical axis
EXTERNC MIA_RESULT_CODE IPL_GEOM_MirrorX(
  unsigned char* dst,           // destination
  unsigned int dststr,         // stride
  const unsigned char* src,     // source
  unsigned int srcstr,         // stride
  unsigned int w,              // ROI size
  unsigned int h);             //

// mirror image relative to vertical axis 
// different BpP values possible
EXTERNC MIA_RESULT_CODE IPL_GEOM_MirrorX2(
  unsigned char* dst,           // destination
  const unsigned char* src,     // source
  int xs,                       // image size
  int ys,                       //
  int BpP);                     // bytes per pixel

EXTERNC void IPL_GEOM_MirrorX3(
  unsigned char* im,        // destination
  unsigned int xs,          // image size
  unsigned int ys);         //

// mirror image relative to horizontal axis 
// different BpP values possible
EXTERNC MIA_RESULT_CODE IPL_GEOM_MirrorY(
  unsigned char* dst,           // destination
  const unsigned char* src,     // source
  int xs,                       // image size
  int ys,                       //
  int BpP);                     // bytes per pixel

// flip image 90 degrees counterclockwise
// BpP!=1 can be used
EXTERNC MIA_RESULT_CODE IPL_GEOM_Flip90(
  unsigned char* dst,           // destination
  const unsigned char* src,     // source
  int xs,                       // image size
  int ys,                       //
  int BpP);                     // bytes per pixel

//== histogram functions ====================================================
EXTERNC MIA_RESULT_CODE IPL_HIST_CalcMedianInLine(
  unsigned char* line,
  const unsigned char* im,
  int xs,
  int ys,
  int h,
  int kew);

// Returns: median
EXTERNC int IPL_HIST_grouping(
  short *bufin,
  int med,
  int lstr);

// Returns: median
EXTERNC int IPL_HIST_grouping1(
  short *bufin,
  int med,
  int lstr);

// make median of all array for byte
EXTERNC int IPL_HIST_mdn_pixE(
  const unsigned char* bufin, // input array
  int                  lstr);  // size bufin

EXTERNC int IPL_HIST_mdn_pixS(
  const short* bufin, // input array
  int          lstr);  // size bufin

// get index of the histogram, below which a given proportion of mass resides
EXTERNC int IPL_HIST_GetQuantile(     // RET: position of quantile
  const unsigned int* pnHist, // IN:  histogram
  unsigned int nMass,         // IN:  histogam mass
  unsigned int nQuantile);    // IN:  quantile in 1/1000

EXTERNC int IPL_HIST_GetQuantileExt(  // RET: position of quantile
  const unsigned int* pnHist, // IN:  histogram
  int nLen,           // IN:  histogram length
  int nMass,          // IN:  histogam mass
  int nQuantile);     // IN:  quantile in 1/1000

// blur histogram
EXTERNC void IPL_HIST_Blur(
  int* dst,         // OUT: blurred hist
  const int* val,   // IN:  source hist
  int len,          // IN:  length of the hist
  int hw);          // IN:  half-window

// blur histogram
EXTERNC void IPL_HIST_Blur_double(
  double* dst,          // OUT: blurred hist
  const double* val,    // IN:  source hist
  int len,              // IN:  length of the hist
  int hw);              // IN:  half-window

// needed for filling brightness histogram gaps
// occurring after color adjustment procedures
// uses linear interpolation by non-zero neighbors
EXTERNC void IPL_HIST_FillGaps(
  int* dst,     // IN/OUT: histogram (transformed in-place)
  int len,      // IN:  histogram data length
  int wid);     // IN:  maximum width of gap

// normalize histogram to a given value
EXTERNC MIA_RESULT_CODE IPL_HIST_Normalize(
  int* dst,         // OUT: normalized histogram
  const int* src,   // IN:  source histogram, if NULL, transform is in-place for dst
  int len,          // IN:  number of elements
  int maxabsval);   // IN:  maximum absolute value

EXTERNC MIA_RESULT_CODE IPL_HIST_EqualizeHist(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys);

EXTERNC MIA_RESULT_CODE IPL_HIST_EqualizeHistExt(
  unsigned char* dst,
  int dststr,
  const unsigned char* src,
  int srcstr,
  int xs,
  int ys);

// 1D dilation, source should have only two values: 
// zero and non-zero (equal for all non-zero points)
EXTERNC void IPL_Morph_DilateLineBinaryCycled(
  int* dst,
  const int* src,
  int len,
  int wid);

// 1D erosion, source should have only two values: 
// zero and non-zero (equal for all non-zero points)
EXTERNC void IPL_Morph_ErodeLineBinaryCycled(
  int* dst,
  const int* src,
  int len,
  int wid);

//== pumps (changing images locally, unlike GEOM) ==============================
// interpolate odd rows with even rows / or vice versa
// BpP!=1 can be used
EXTERNC MIA_RESULT_CODE IPL_PMP_Interpolate(
  unsigned char* dst,           // destination
  const unsigned char* src,     // source
  int xs,                       // image size
  int ys,                       //
  int BpP,                      // bytes per pixel
  int bInterpolateOdd);         // odd = f(even)

// copy odd rows with even rows / or vice versa
// BpP!=1 can be used
EXTERNC MIA_RESULT_CODE IPL_PMP_CopyImageStrings(
  unsigned char* dst,           // destination
  const unsigned char* src,     // source
  int xs,                       // image size
  int ys,                       //
  int BpP,                      // bytes per pixel
  int bInterpolateOdd);         // odd = f(even)

EXTERNC void IPL_PMP_PyramidMidOf8x4IL(
  unsigned char* dst,
  const unsigned char* src,
  int xs,   // source size
  int ys);

EXTERNC void IPL_PMP_PyramidMidOf4x2IL(
  unsigned char* dst,
  const unsigned char* src,
  int xs,   // source size
  int ys);

// make transformation of 4x4 square to one pixel
// even interlace half-frame is taken, so 4*2 points are actually used
EXTERNC void IPL_PMP_PyramidMidOf4x2IL_int(
  int* dst,
  const unsigned char* src,
  int xs,   // source size
  int ys);

EXTERNC void IPL_PMP_PyramidMidOf2x1IL(
  unsigned char* dst,
  const unsigned char* src,
  int xs,   // source size
  int ys);

EXTERNC void IPL_PMP_PyramidMidOf2x2(
  unsigned char* dst,
  const unsigned char* src,
  int xs,   // source size
  int ys);

// summing downscaling by 2
EXTERNC MIA_RESULT_CODE IPL_PMP_Squeeze2_Good(
  unsigned char* dst,
  int dststr,
  const unsigned char* src,
  int srcstr,
  int xs,
  int ys,
  int BpP);

// simple downscaling by 2
EXTERNC MIA_RESULT_CODE IPL_PMP_Squeeze2_Simple(
  unsigned char* dst,
  int dststr,
  const unsigned char* src,
  int srcstr,
  int xs,
  int ys,
  int BpP);

// interpolating upscaling by 2
EXTERNC MIA_RESULT_CODE IPL_PMP_Expand2_Good(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  int BpP);

// simple upscaling by 2
EXTERNC MIA_RESULT_CODE IPL_PMP_Expand2_Simple(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  int BpP);

// combining half-frames
EXTERNC MIA_RESULT_CODE IPL_PMP_CombineHalfFrames(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  int bOddIsFirst);

// splitting half-frames
EXTERNC MIA_RESULT_CODE IPL_PMP_SplitHalfFrames(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  int bOddIsFirst);

// crop piece of source image to destination image
EXTERNC void IPL_PMP_Crop(
  unsigned char* dst,         // IN:  destination image
  int dstr,                   // IN:  dst image stride
  int dxs,                    // IN:  dst image size
  int dys,                    //
  int xplace,                 // IN:  position to insert the ROI
  int yplace,                 //
  const unsigned char* src,   // IN:  source image
  int sstr,                   // IN:  src image stride
  int sxs,                    // IN:  src image size
  int sys,                    //
  int x,                      // IN:  source ROI
  int y,                      //
  int w,                      //
  int h);                     //

// crop piece of source image to destination image
EXTERNC void IPL_PMP_CropCentering(
  unsigned char* dst,         // IN:  destination image
  int dxs,                    // IN:  dst image size
  int dys,                    //
  const unsigned char* src,   // IN:  source image
  int sxs,                    // IN:  src image size
  int sys,                    //
  int xc,                     // IN:  center coordinates
  int yc);                    //

EXTERNC void IPL_PMP_MakeIntegralImage(
  int* intim,               // OUT: integer image
  const unsigned char* im,  // IN:  source image
  int xs,                   // IN:  size
  int ys);

EXTERNC void IPL_PMP_MakeIntegralImage_prim(
  int* intim,               // OUT: integral image, stride is (xs+1) elements
  const unsigned char* im,  // IN:  source image
  int xs,                   // IN:  size
  int ys);

EXTERNC MIA_RESULT_CODE IPL_PMP_Normalize(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  int lo_cut,     // low quantile cutoff, in 1/1000 shares
  int hi_cut);    // high quantile cutoff, in 1/1000 shares

EXTERNC MIA_RESULT_CODE IPL_PMP_ScaleYUVtoGRAY(
  unsigned char* dst, // destination image
  int dststr,         // stride in bytes
  int W,              // destination image size
  int H,              //
  const unsigned char* src,   // source image
  int srcstr,         // stride in bytes
  int xs,             // source image size
  int ys);            //

EXTERNC MIA_RESULT_CODE IPL_PMP_ScaleBAYERtoGRAY(
  unsigned char* dst, // destination image
  int dststr,         // stride in bytes
  int W,              // destination image size
  int H,              //
  const unsigned char* src,   // source image
  int srcstr,         // stride in bytes
  int xs,             // source image size
  int ys);             //

EXTERNC MIA_RESULT_CODE IPL_PMP_ScaleBAYERtoGRAYcrop(
  unsigned char* dst, // destination image
  int dststr,         // stride in bytes
  int W,              // destination image size
  int H,              //
  const unsigned char* src,   // source image
  int srcstr,         // stride in bytes
  int xs,             // source image size
  int ys,
  int roix,                   // source ROI : x
  int roiy,					  //              y
  int roiw,                   // ROI width
  int roih);    			  //	 height

EXTERNC MIA_RESULT_CODE IPL_PMP_ScaleRGBtoGRAY(
  unsigned char* dst, // destination image
  int dststr,         // stride in bytes
  int W,              // destination image size
  int H,              //
  const unsigned char* src,   // source image
  int srcstr,         // stride in bytes
  int xs,             // source image size
  int ys);            //

EXTERNC MIA_RESULT_CODE IPL_PMP_ScaleGrayToGray(
  unsigned char* dst, // destination image
  int dststr,         // stride in bytes
  int W,              // destination image size
  int H,              //
  const unsigned char* src,   // source image
  int srcstr,         // stride in bytes
  int xs,             // source image size
  int ys);            //

// make time resampling transformation
EXTERNC MIA_RESULT_CODE IPL_PMP_ResampleImageSequence(
  unsigned char* pDst,        // OUT: resampled sequence
  int nFrDst,                 // IN:  number of full frames in the resampled sequence
  const char* pBrtShifts,     // IN:  brightness shifts, may be NULL
  const unsigned char* pSrc,  // IN:  source image sequence
  int nFrSrc,                 // IN:  number of full frames in source
  int nInserted,              // IN:  number of duplicated half-frames in the beg.
  int xs,                     // IN:  image width in pixels
  int ys,                     // IN:  image height in pixels
  int bOddIsFirst);           // IN:  order of half-frames

EXTERNC void IPL_PMP_SmoothX_short(
  short* sm1,
  const unsigned char* im,
  int xs,
  int ys);

EXTERNC void IPL_PMP_SmoothY_short(
  short* sm1,
  const short* im,
  int xs,
  int ys);

EXTERNC void IPL_PMP_SmoothX_int(
  int* sm1,
  const int* im,
  int xs,
  int ys);

EXTERNC void IPL_PMP_SmoothX_intchar(
  int* sm1,
  const uint8* im,
  int xs,
  int ys);

EXTERNC void IPL_PMP_SmoothY_int(
  int* sm1,
  const int* im,
  int xs,
  int ys);

//== projections =========================================================
// calculate left side, right side and total circular projection 
// in a concentric ring
EXTERNC MIA_RESULT_CODE IPL_PROJ_FindHoughDonatorProjection3(
  int* pnProjT,             // OUT: circular projection - total
  int* pnProjR,             // OUT: circular projection - right side
  int* pnProjL,             // OUT: circular projection - left side
  const unsigned char* im,  // IN:  image
  int xs,                   // IN:  image size
  int ys,                   // 
  int xc,                   // IN:  ring center
  int yc,                   //
  int r,                    // IN:  inner radius
  int* pR,                  // IN/OUT: outer radius proposed/actually used
  int nMaxRadOfDence,       // IN:  maximum radius for dence processing. 0 - dence always
  void* buf,                // IN:  external buffer
  int* buflen);             // IN/OUT: allocated/used bytes

// calculate four quadrant projections (for pupil-from-iris)
EXTERNC MIA_RESULT_CODE IPL_PROJ_FindHoughDonatorProjection4(
  int* pnProjR,             // OUT: circular projection - right side
  int* pnProjL,             // OUT: circular projection - left side
  int* pnProjT,             // OUT: circular projection - top side
  int* pnProjB,             // OUT: circular projection - bottom side
  const unsigned char* im,  // IN:  image
  int xs,                   // IN:  image size
  int ys,                   // 
  int xc,                   // IN:  center
  int yc,                   //
  int R,                    // IN:  outer radius
  void* buf,                // IN:  temporary buffer
  int* plen,                // IN/OUT: temporary buffer allocated/used
  const char* nambeg);

// calculate left and right side projection histograms in a concentric ring
EXTERNC MIA_RESULT_CODE IPL_PROJ_FindHoughDonatorProjection5(
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
  int* buflen);             // IN/OUT: allocated/used bytes

// same as '5', but with thresholds from '3'
MIA_RESULT_CODE IPL_PROJ_FindHoughDonatorProjection6(
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
  int* buflen);             // IN/OUT: allocated/used bytes

// same as '5' but with blurring
MIA_RESULT_CODE IPL_PROJ_FindHoughDonatorProjection7(
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
  int* buflen);             // IN/OUT: allocated/used bytes

//== pyramidal functions ==================================================
// estimate image sharpness by pyramidal fractal dimension
EXTERNC MIA_RESULT_CODE IPL_PYR_EstimateSharpness(
  int* pnSharpness,         // OUT: sharpness estimation
  const unsigned char* im,  // IN:  source image
  int xs,                   // IN:  source image size
  int ys,                   //
  void* pvTmpBuf,           // IN:  external buffer
  int* pnBufSiz);           // IN/OUT: bytes alloced/used in buffer

// estimate image sharpness by pyramidal fractal dimension
EXTERNC MIA_RESULT_CODE IPL_PYR_EstimateSharpness_withMask(
  int* pnSharpness,           // OUT: motion estimation
  const unsigned char* im,    // IN:  source image
  const unsigned char* mask,  // IN:  mask
  int xs,                     // IN:  source image size
  int ys,                     //
  void* pvTmpBuf,             // IN:  external buffer
  int* pnBufSiz);             // IN/OUT: bytes alloced/used in buffer

// estimate image motion as difference of half-frames
EXTERNC MIA_RESULT_CODE IPL_PYR_EstimateMotion(
  int* pnMotion,            // OUT: motion estimation
  const unsigned char* im,  // IN:  source image
  int xs,                   // IN:  source image size
  int ys,                   //
  void* pvTmpBuf,           // IN:  external buffer
  int* pnBufSiz);           // IN/OUT: bytes alloced/used in buffer

// estimate image noise as a difference between original and its median
MIA_RESULT_CODE IPL_PYR_EstimateNoise(
  int* pnNoise,             // OUT: noise estimation
  const unsigned char* im,  // IN:  source image
  unsigned int xs,                   // IN:  source image size
  unsigned int ys,                   //
  void* pvTmpBuf,           // IN:  external buffer
  int* pnBufSiz);           // IN/OUT: bytes alloced/used in buffer

// detect loss of sync as a spike of differences between adjacent lines
EXTERNC MIA_RESULT_CODE IPL_PYR_DetectSyncLoss(
  int* pnSyncLoss,          // OUT: sync loss detection flag
  const unsigned char* im,  // IN:  source image
  int xs,                   // IN:  source image size
  int ys,                   //
  void* pvTmpBuf,           // IN:  external buffer
  int* pnBufSiz);           // IN/OUT: bytes alloced/used in buffer

//== image segmentation ======================================================
EXTERNC MIA_RESULT_CODE IPL_SEGM_FindMaximaDigging_uint8(
  SValuePosition* maxima,       // OUT:    destination array
  int* size,           // IN/OUT: size of array in elements/number of extremums found
  uint8* array,         // IN:     image to be searched for maxima (zeroed on exit)
  int w,                        // IN:     image size
  int h,                        //
  uint8 maxthreshold,   // IN:     threshold of max
  void* extbuftmp);             // IN:     temporary buffer

EXTERNC MIA_RESULT_CODE IPL_SEGM_FindMaximaDigging_uint8_ext(
  SValuePosition* maxima,   // OUT:    destination array
  int* size,        // IN/OUT: size of array in elements/number of extremums found
  uint8* array,     // IN:     image to be searched for maxima (zeroed on exit)
  int w,            // IN:     image size
  int h,            //
  uint8 minval_thr, // IN: minimum possible brightness of local max to be detected
  uint8 mindrop_thr,// IN: minimum possible brightness drop to be detected
  int maxobjsize,   // IN: maximum allowed object size
  void* extbuftmp); // IN:     temporary buffer

EXTERNC MIA_RESULT_CODE IPL_SEGM_FindMaximaDigging_uint32(
  SValuePosition* maxima,       // OUT:    destination array
  int* size,           // IN/OUT: size of array in elements/number of extremums found
  uint32* array,         // IN:     image to be searched for maxima (zeroed on exit)
  int w,                        // IN:     image size
  int h,                        //
  uint32 maxthreshold,   // IN:     threshold of max
  void* extbuftmp);             // IN:     temporary buffer

// clusterize elements with minimal-distance joining rule
EXTERNC MIA_RESULT_CODE IPL_SEGM_ClusteriseByMinDist(
  int* anClusterNumbers,    // OUT: cluster numbers
  int* pnClasses,           // IN/OUT: number of clusters required or
                            //         zero to judge itself/obtained
  int nElements,            // IN:  number of elements
  const void** pavFeatures, // IN:  features to be compared
  unsigned int (*pfnCompareFeatures)(const void*,const void*), 
                            // IN:  comparator. Must return -1 for infinity
  unsigned char* pucBuffer, // IN:  external memory buffer
  int* pnBufsize);          // IN/OUT: buffer size allocated/used

EXTERNC MIA_RESULT_CODE IPL_SEGM_KillBlik(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  void* buf,
  int* len);

//== Viola-Jones cascade ======================================================
typedef struct
{
  double weight;  // weight
  int nType;      // type of primitive
  int nThreshold; // classifier threshold
  int nThreshold2; // classifier threshold
//  int x[4];       // positions of window parts (not all used)
  //int y[4];       //
  int idx[8];
} SViolaJonesPrimitive;

typedef struct
{
  int nPrimitives;  // number of primitives
//  double* pWeights; // weights of primitives
  SViolaJonesPrimitive* pPrimitives;
} SViolaJonesCascade;

typedef struct
{
  int xs; // size of window
  int ys; // 
  int nCascade; // number of cascades
  SViolaJonesCascade* pCascades;
} SViolaJonesClassifier;


EXTERNC MIA_RESULT_CODE VJ_train(
  SViolaJonesClassifier** ppSVJ,
  const uint8* pPosSamples,
  int nPosSamples,
  const uint8* pNegSamples,
  int nNegSamples,
  int xs,
  int ys,
  const char* repa,
  const char* nambeg);

EXTERNC MIA_RESULT_CODE VJ_Save(
  const SViolaJonesClassifier* pS,
  const char* nama);

EXTERNC MIA_RESULT_CODE VJ_Load(
  const SViolaJonesClassifier** ppS,
  const char* nama);

EXTERNC MIA_RESULT_CODE VJ_Free(
  SViolaJonesClassifier** ppS);

EXTERNC MIA_RESULT_CODE VJ_Classify(
  int* pDecision,
  const SViolaJonesClassifier* pSVJ,
  const uint8* im,
  int str);

//== CVL remnants =============================================================
// expands rect ROI on 2*HalfWidthSizeExpand and 2*HalfHeightSizeExpand values
EXTERNC void CVL_ExpandRectROI(
  int*                pDestROIStartX,        // Dest ROI left
  int*                pDestROIStartY,        // Dest ROI top
  unsigned int*       pDestROIWidth,         // Dest ROI width  
  unsigned int*       pDestROIHeight,        // Dest ROI height
  int                 BoundROIStartX,        // Bounding ROI left
  int                 BoundROIStartY,        // Bounding ROI top
  unsigned int        BoundROIWidth,         // Bounding ROI width 
  unsigned int        BoundROIHeight,        // Bounding ROI height
  const unsigned int  HalfWidthSizeExpand,   // Half of width expand value
  const unsigned int  HalfHeightSizeExpand); // Half of height expand value

EXTERNC void IPL_POLAR_TransformBilinearHorFillByLast(
  unsigned char*       Dest, 
  const unsigned char* Source,
  unsigned int         DestWidth,
  unsigned int         DestHeight,
  unsigned int         DestWidthItems,
  unsigned int         SourceWidth,
  unsigned int         SourceHeight,
  unsigned int         SourceWidthItems,
  int                  SrcCx, 
  int                  SrcCy,
  unsigned int         SrcStartRadius,
  unsigned int         SrcFinishRadius,
  float                SrcStartAngleRad,
  float                SrcFinishAngleRad,
  const unsigned char  FirstFillValue,
  int*                 HorizontalBuffer,    // SizeItems=DestWidth
  float*               VerticalBuffer);     // SizeItems=DestHeight

EXTERNC void IPL_POLAR_MakePolarRubberInt(
  uint8* polar,
  int32 w,
  int32 numLines,
  const uint8* im,
  int32 xs,
  int32 ys,
  int32 xC, 
  int32 yC, 
  int32 rC, 
  int32 xCI, 
  int32 yCI, 
  int32 rCI);

// make polar rubber sheet transform
EXTERNC void IPL_POLAR_MakePolarRubber(
  uint8* polar,
  int32 w,
  int32 numLines,
  const uint8* im,
  int32 xs,
  int32 ys,
  int32 xC, 
  int32 yC, 
  int32 rC, 
  int32 xCI, 
  int32 yCI, 
  int32 rCI);

EXTERNC void MakePolarAllINCITS2(
  uint8* polar,
  int32 w,
  int32 numLines,
  const uint8* im,
  int32 xs,
  int32 ys,
  int32 xC, 
  int32 yC, 
  int32 rC, 
  int32 xCI, 
  int32 yCI, 
  int32 rCI);     // SizeItems=DestHeight

EXTERNC void IPL_POLAR_MakePolarRubberDeform(
  uint8* polar,
  int32 w,
  int32 numLines,
  const uint8* im,
  int32 xs,
  int32 ys,
  int32 xC, 
  int32 yC, 
  int32 rC, 
  int32 xCI, 
  int32 yCI, 
  int32 rCI,
  int32 DeformMethod, // deformation method; 0-usual rubber; 1-power mode; 2-ivan deformation 
  double dscaler);     // parameter of deformation for method 1 and 2

EXTERNC void CVL_EdgeSobelFast_ByteToUShort(                          
  unsigned short*       Dest,          // destination image
  const unsigned char*  Source,        // source image
  const unsigned int    Width,         // width of images (in elements)
  const unsigned int    Height,        // height of images (in elements)
  const unsigned int    Stride,        // row stride of image (in elements)
  unsigned int          ROIStartX,     // x position of ROI (in elements)
  unsigned int          ROIStartY,     // y position of ROI (in elements)
  unsigned int          ROIWidth,      // width of ROI (in elements)
  unsigned int          ROIHeight,     // height of ROI (in elements)
  int*                  MemoryBufferX, // SizeItems = 3*ROIWidth;
  int*                  MemoryBufferY);// SizeItems = 3*ROIWidth;

EXTERNC void CVL_EdgeSobelFast_Byte( 
  unsigned char*        Dest,           // destination image
  const unsigned char*  Source,         // source image
  const unsigned int    Width,          // width of images (in elements)
  const unsigned int    Height,         // height of images (in elements)
  const unsigned int    Stride,         // row stride of image (in elements)
  unsigned int          ROIStartX,      // x position of ROI (in elements)
  unsigned int          ROIStartY,      // y position of ROI (in elements)
  unsigned int          ROIWidth,       // width of ROI (in elements)
  unsigned int          ROIHeight,      // height of ROI (in elements)
  int*                  MemoryBufferX,  // SizeItems = 3*ROIWidth;
  int*                  MemoryBufferY); // SizeItems = 3*ROIWidth;

EXTERNC void CVL_DownScaleImage2TimesByCalcAverageInRect(
  unsigned char*       Dest,              //Destination
  const unsigned char* Source,            //Source
  unsigned int         DestROIStartX,
  unsigned int         DestROIStartY,
  unsigned int         SourceROIStartX,   // Source ROI
  unsigned int         SourceROIStartY,   //
  unsigned int         SourceROIWidth,    //
  unsigned int         SourceROIHeight,   //
  unsigned int         DestWidthBytes,    //Row stride of destination image
  unsigned int         SourceWidthBytes); //Row stride of source image

#define ARRAY_ITEM_ALIGN_POWER (4)
#define ALIGN_IMAGE_POWER ARRAY_ITEM_ALIGN_POWER
#define VALUE_IMAGE_TO_ADD   ((1<<(ALIGN_IMAGE_POWER)) - 1)
#define VALUE_IMAGE_TO_AND   (~VALUE_IMAGE_TO_ADD)

#define GET_IMAGE_WIDTH_ITEMS(nWidth,Channels)\
  ((unsigned int)((nWidth * Channels + VALUE_IMAGE_TO_ADD) & VALUE_IMAGE_TO_AND))

typedef struct
{
  unsigned char* Pixels;
  unsigned int   Width;
  unsigned int   Height;
  unsigned int   WidthItems;
  unsigned int   Channels;
} CVL_IMAGE;

typedef struct
{
  unsigned short* Pixels;
  unsigned int    Width;
  unsigned int    Height;
  unsigned int    WidthItems;
  unsigned int    Channels;
} CVL_ushIMAGE;

EXTERNC void CVL_BilinearStretchInROIFill(
  unsigned char*       Dest,                // Destination
  const unsigned char* Source,              // Source
  unsigned int         DestWidth,           // Width of Dest image.
  unsigned int         DestHeight,          // Height of Dest image.
  unsigned int         DestWidthBytes,      // Dest image stride
  unsigned int         SourceWidth,         // Width of Source image.
  unsigned int         SourceHeight,        // Height of Source image.
  unsigned int         SourceWidthBytes,    // Source image stride
  int                  SourceROIStartX,     // ROI on SOURCE image. On DEST image it is 0.
  int                  SourceROIStartY,     // ROI on SOURCE image. On DEST image it is 0.
  unsigned int         SourceROIWidth,      // ROI on SOURCE image. On DEST image it is DestWidth.
  unsigned int         SourceROIHeight,     // ROI on SOURCE image. On DEST image it is DestHeight.
  unsigned int*        HorizontalBuffer,    // SizeItems=DestWidth
  unsigned int*        VerticalBuffer);     // SizeItems=DestHeight

EXTERNC void CVL_FillCircleInRectUShort(
  unsigned short* Source,
  unsigned int   Width,
  unsigned int   Height,
  unsigned int   WidthItems,
  unsigned int   ROIStartX,
  unsigned int   ROIStartY,
  unsigned int   ROIWidth,
  unsigned int   ROIHeight,
  int            XCenter,
  int            YCenter,
  unsigned int   Radius,
  unsigned short ValueToSet);

EXTERNC void CVL_FillHorHalfCircleRingInRectUShort(
  unsigned short* Source,
  unsigned int   Width,
  unsigned int   Height,
  unsigned int   WidthItems,
  unsigned int   ROIStartX,
  unsigned int   ROIStartY,
  unsigned int   ROIWidth,
  unsigned int   ROIHeight,
  int            XCenter,
  int            YCenter,
  unsigned int   StartRadius,
  unsigned int   FinalRadius,
  unsigned short ValueToSet);

EXTERNC void CVL_CalcMaxGradAreaUShort(
  unsigned int*         pMaxGradArea,
  unsigned int*         pGradArea,
  const unsigned short* Source,                   
  unsigned int          ROIStartX,
  unsigned int          ROIStartY,
  unsigned int          ROIWidth,
  unsigned int          ROIHeight,
  unsigned int          WidthItems);

EXTERNC void CVL_CalcMaxGradAreaInCircHorHalfRingUShort(
  unsigned int*         pMaxGradArea,
  unsigned int*         pGradArea,
  const unsigned short* Source,                   
  unsigned int          ROIStartX,
  unsigned int          ROIStartY,
  unsigned int          ROIWidth,
  unsigned int          ROIHeight,
  unsigned int          WidthItems,
  int                   XCenter,
  int                   YCenter,
  unsigned int          StartRadius,
  unsigned int          FinalRadius);

typedef struct
{
  int nS;
  int nI;
} SBlobData1;

typedef struct
{
  int nParent;
  int nChild1;
  int nChild2;
  int nNeibStart;
  int nNeibCount;
} SBlobInfo;

typedef struct
{
  int64 nS;
  int64 nSx;
  int64 nSy;
  int64 nSxx;
  int64 nSxy;
  int64 nSyy;
  int64 nSgx;
  int64 nSgy;
} SBlobData3;

typedef struct
{
  int64 nS;
  int64 nSx;
  int64 nSy;
  int64 nSxx;
  int64 nSxy;
  int64 nSyy;
  int64 nSxxx;
  int64 nSxxy;
  int64 nSxyy;
  int64 nSyyy;
  int64 nSxxxx;
  int64 nSxxyy;
  int64 nSyyyy;
} SBlobData4;

typedef struct 
{
  SBlobData4 moments;
  double x0;  // approximation
  double y0;
  double r;
  double q;
  int  nPointCount; // point list
  int* pnPointList;
  int  nNeibCount;  // neighbors list
  int* pnNeibList;
} SArchunkData2;

EXTERNC MIA_RESULT_CODE IPL_ELL_GetApproximatingCircle(
  double *px0,
  double *py0,
  double *pr,
  const SBlobData4* pD);

// up-to-low recursion
EXTERNC MIA_RESULT_CODE IPL_CONOB_GetBlobPointIndices(
  int** ppnPointList,       // OUT: list of extracted points
  int** ppnPointCounts,     // OUT: triples of <start_position,length,src_blob>
  int*  pnBlobCount,        // OUT: number of resulting blobs
  const SBlobInfo *psBI,    // IN:  source blob infos
  int nBlobSrc);            // IN:  number of source blobs (in blob heap)

// collect points to handy array
EXTERNC MIA_RESULT_CODE IPL_CONOB_GetBlobPointCoords(
  int** ppnPointCoords,     // OUT: coordinates as pairs <x,y>
  const SBlobData4* psBD,   // IN:  blob data
  const int* pnPointList,   // IN:  list of point indexes
  const int* pnPointCounts, // IN:  triples of <start_position,length,src_blob>
  int nPoints);

EXTERNC MIA_RESULT_CODE IPL_CONOB_Enumerate_v1(
  SBlobInfo **ppsBI,        // OUT: blob infos
  SBlobData1 **ppsBD,       // OUT: blob data
  int **ppnNeibList,        // OUT: neighbor indices
  int *pnBlob,              // OUT: number of blobs
  const unsigned char* im,  // IN:  source image
  int xs,                   // IN:  image size
  int ys,                   //
  const char* nambeg);      // IN:  for debug, set to NULL

EXTERNC MIA_RESULT_CODE IPL_CONOB_Enumerate_v2(
  SBlobInfo **ppsBI,        // OUT: blob infos
  SBlobData1 **ppsBD,       // OUT: blob data
  int **ppnNeibList,        // OUT: neighbor indices
  int *pnBlob,              // OUT: number of blobs
  const unsigned char* im,  // IN:  source image
  int xs,                   // IN:  image size
  int ys,                   //
  const char* nambeg);      // IN:  for debug, set to NULL

EXTERNC MIA_RESULT_CODE IPL_CONOB_Enumerate_StraightChunks_v1(
  SBlobInfo **ppsBI,        // OUT: blob infos
  SBlobData3 **ppsBD,       // OUT: blob data
  int **ppnNeibList,        // OUT: neighbor indices
  int *pnBlob,              // OUT: number of blobs
  const unsigned char* im,  // IN:  source image
  const short* pgx,
  const short* pgy,
  int xs,                   // IN:  image size
  int ys,                   //
  const char* nambeg);      // IN:  for debug, set to NULL

EXTERNC MIA_RESULT_CODE IPL_CONOB_Enumerate_CircularChunks_v1(
  SBlobInfo  **ppsBI,       // OUT: blob infos
  SBlobData4 **ppsBD,       // OUT: blob data
  int **ppnNeibList,        // OUT: neighbor indices
  int *pnBlob,              // OUT: number of blobs
  const unsigned char* im,  // IN:  source image
  const short* pgx,
  const short* pgy,
  int xs,                   // IN:  image size
  int ys,                   //
  const char* nambeg);      // IN:  for debug, set to NULL

#endif /* __ipl_h__ */
