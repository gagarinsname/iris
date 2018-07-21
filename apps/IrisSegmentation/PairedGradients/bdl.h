/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  23 April 2002                                          */
/*      Modified: many times                                             */
/*      Revision: x.0.00                                                 */
/*      Purpose:  Iris library API (common part)                         */
/*      Authors:                                                         */
/*        Evgeny Demin                                                   */
/*        Konstantin Gankin                                              */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#ifndef __bdl_h__
#define __bdl_h__

#include "errorcodes.h"
#include "stddefs.h"
#include "ipl.h"

EXTERNC MIA_RESULT_CODE BDL_GetVersionInfo(
  char* pcVersionString,  // OUT: place reserved for asciiz string
  int pnStringSize);      // IN:  string length

//== data types and structures ==

//== center detection types
#define IRL_CENTER_FLAG_CALC_FAST      0x00000001 // fast error-prone calculation
#define IRL_CENTER_FLAG_CALC_ROBUST    0x00000002 // slow reliable calculation
#define IRL_CENTER_FLAG_CALC_BIPOINT   0x00000004 // bi-point method
#define IRL_CENTER_FLAG_CALC_CANNY     0x00000010 // use Canny pre-selection
#define IRL_CENTER_FLAG_CALC_BLOB      0x00000020 // blob method
#define IRL_CENTER_FLAG_CALC_RELAX     0x00000040 // use relaxation
#define IRL_CENTER_FLAG_CALC_FLASHMASK 0x00000080 // use flash mask
#define IRL_CENTER_FLAG_CALC_MSD       0x00000100 // use MSD method
//== center detection options
#define IRL_CENTER_FLAG_EVAL_QUALITY   0x00010000 // evaluate quality (for robust)
#define IRL_CENTER_FLAG_NORM_BRT       0x00020000 // normalize brightness (for robust)
#define IRL_CENTER_FLAG_NORM_GRAD      0x00040000 // normalize gradient (for robust)
#define IRL_CENTER_FLAG_ERODE_VOTE     0x00080000 // erode color voting value (for robust)
#define IRL_CENTER_FLAG_ROI20          0x00100000 // search center in 20% region
#define IRL_CENTER_FLAG_TINYPUPIL      0x00200000 // very small pupil
#define IRL_CENTER_FLAG_SCALE_FACTOR   0x10000000 // two highest bits set scale factor (for robust)
//== projection-relevant functions options
#define IVIR_ARRPOXIR_LEFEYE  0x00000001  // definitely left eye
#define IVIR_ARRPOXIR_RIGEYE  0x00000002  // definitely right eye
#define IVIR_ARRPOXIR_GENPROJ 0x00100000  // generate projections
#define IVIR_ARRPOXIR_USEPROJ 0x00200000  // use given projections
#define IVIR_ARRPOXIR_ADDHEUR 0x00010000  // additional heuristics
#define IVIR_ARRPOXIR_PRPRPR  0x00020000  // PRePRocess PRojection

typedef struct
{
  int xc;   // X coordinate in fix20 format
  int yc;   // Y coordinate in fix20 format
  int q;    // quality in scale 0-1000
} SCenterInfo;

typedef struct
{
  int pu_xc;    // X coordinate of supposed pupil
  int pu_yc;    // Y coordinate of supposed pupil
  int pu_r;     // supposed pupil radius
  int pu_q;     // quality in -1 - no attempt to detect, 0 - not detected, 1-100
  int ir_xc;    // X coordinate of supposed iris
  int ir_yc;    // Y coordinate of supposed iris
  int ir_r;     // supposed iris radius (BazRad)
  int ir_q;     // quality in -1 - no attempt to detect, 0 - not detected, 1-100
} SBazRadInfo;

typedef struct
{
  int xc;   // X coordinate
  int yc;   // Y coordinate
  int r;    // radius (elliptic big axis)
  int q1;   // quality?
  int q2;   // quality?
} SPupilInfo;

typedef struct
{
  int xc;   // X coordinate
  int yc;   // Y coordinate
  int r;    // radius
  int q1;   // quality?
  int q2;   // quality?
  int OccUpB; // beginning of upper occlusion in grads
  int OccUpE; // ending of upper occlusion in grads
              // counted counterclockwise from 3hrs
  int OccDnB; // beginning of lower occlusion in grads
  int OccDnE; // ending of lower occlusion in grads
              // counted clockwise from 3hrs
} SIrisInfo;

typedef struct
{
  int x;          // position of center
  int y;
  int a;          // major axis
  int b;          // minor axis
  int t;          // tilt in FIX20
  int qLocation;  // location quality
  int qGoodness;  // goodness quality
} SEllipseData;

typedef struct
{
  int fl1;          // is flash 1 off/on (0/1)
  int fl2;          // is flash 2 off/on (0/1)
  int fl3;          // is flash 3 off/on (0/1)
  int fl4;          // is flash 4 off/on (0/1)
} SFlashOnInfo;

typedef struct
{
  SEllipseData PupEll;
  SEllipseData IrEll;
} SEyeEllipseData;

#define DET_BORD_RADS 360

typedef struct
{
  int          xc;
  int          yc;
  unsigned int Count;     // number of border points
  int          Radii[DET_BORD_RADS];     // radii
  int          Qualities[DET_BORD_RADS]; // qualities
  int          bValid;    // overall detection quality
} SPolarPointData;

typedef struct
{
  int            UpperAngleBeg;
  int            UpperAngleEnd;
  int            LowerAngleBeg;
  int            LowerAngleEnd;
} SOcclusionData;

typedef struct
{
  SPolarPointData  pupil;
  SPolarPointData  iris;
  SOcclusionData   occlusion;
  int              bValid;              // overall detection quality
} SEyeRefinedData;

typedef struct
{
  // source data
  unsigned char* im;  // image, its size and stride
  int xs;
  int ys;
  int bInterlaced;
  // intermediate data
  SCenterInfo sCI;
  SBazRadInfo sBRI;
  SPupilInfo  sPI;
  SIrisInfo   sII;
  SEyeEllipseData eye_ellipse;
  SEyeRefinedData ERD;
  unsigned char* TotalMask;
  unsigned char* FlashMask;
  unsigned char* OcclusionMask;
  int nTotalMaskQuality;  // percent of coverage of iris with mask
  int nDetQuality;
  // debug data
  char* acFileNameBeginning;
} SBorderData_ALL;

typedef struct
{ // IQCE qualities
  int nScalarOverallQuality;  // IQCE:Q1
  int nIrisScleraContrast;    // IQCE:Q2
  int nIrisPupilContrast;     // IQCE:Q3
  int nGrayScaleDensity;      // IQCE:Q4
  int nHeadOrientation;       // this is a dummy
  int nSightDirection;        // IQCE:Q6
  int nIrisBoundaryShape;     // IQCE:Q7
  int nIrisSize;              // IQCE:Q8
  int nMotionBlur;            // IQCE:Q9
  int nPupilIrisRatio;        // IQCE:Q10
  int nSharpness;             // IQCE:Q11
  int nSignalToNoiseRatio;    // IQCE:Q12
  int nVisibleIrisArea;       // IQCE:Q13
  int nOcclusions;            // IQCE:Q14
  int nReserved[2];
  int nVendorDefined[16];
} SQualityMeasuresIQCE_v1;    // variant of 'conops_v2'

typedef struct
{ // IQCE qualities
  int nScalarOverallQuality;  // IQCE:Q1
  int nGrayLevelSpread;       // IQCE:Q2  Former Q4:nGrayScaleDensity
  int nIrisSize;              // IQCE:Q3  Former Q8
  int nPupilIrisRatio;        // IQCE:Q4  Former Q10
  int nUsableIrisArea;        // IQCE:Q5  Former Q13:nVisibleIrisArea
  int nIrisScleraContrast;    // IQCE:Q6  Former Q2
  int nIrisPupilContrast;     // IQCE:Q7  Former Q3
  int nIrisShape;             // IQCE:Q8 =new
  int nPupilShape;            // IQCE:Q9  Former Q7:nIrisBoundaryShape
  int nMargin;                // IQCE:Q10=new
  int nSharpness;             // IQCE:Q11
  int nMotionBlur;            // IQCE:Q12 Former Q9
  int nSignalToNoiseRatio;    // IQCE:Q13 Former Q12
  int nMagnification;         // What is it?! No explanation given.
  int nHeadRotation;          // IQCE:Q15 Former nHeadOrientation, calculated now
  int nGazeAngle;             // IQCE:Q16 Former Q6:nSightDirection
  int nInterlace;             // IQCE:Q17=new
  int nReserved[15];
  int nVendorDefined[32];
} SQualityMeasuresIQCE;    // variant of 'conops_v3'

// qualities of individual image defined in ISO 5th working draft
typedef struct
{
  double nFrontalGazeAzimuth;          // Q1
  double nFrontalGazeElevation;        // Q2
  double nGrayScaleUtilization;        // Q3
  double nIrisImageAuthenticity;       // Q4
  double nIrisImageNoise;              // Q5
  double nIrisBoundaryShape;           // Q6
  double nIrisPupilBoundaryContrast;   // Q7
  double nIrisPupilConcentricity;      // Q8
  double nIrisScleraBoundaryContrast;  // Q9
  double nIrisRadius;                  // Q10
  double nMargin;                      // Q11
  double nMotionBlur;                  // Q12
  double nPupilBoundaryShape;          // Q13
  double nPupilIrisRatio;              // Q14
  double nSharpness;                   // Q15
  double nUsableIrisArea;              // Q16
  double nReserved[16];
  double nOverall;             // V1: overall quality
  double nVendorDefined[31];
} SQualityMeasuresISO5;

// parameters of pupil
typedef struct
{
  SPupilInfo pupCircle;
  SEllipseData pupEllipse;
  int det_board_center_x;
  int det_board_center_y;
  short det_bord[DET_BORD_RADS];
  short det_bord_q[DET_BORD_RADS];
  int qUp;  // quality of upper part
  int qLow; // etc...
  int qLeft;
  int qRight;
  int qUpQuarter;
  int qLowQuarter;
  int qLeftQuarter;
  int qRightQuarter;
} BDL_PupilParams;

// parameters of iris
typedef struct
{
  SIrisInfo irisCircle;
  short occl_bord[DET_BORD_RADS];
  short occl_bord_q[DET_BORD_RADS];
  int BangleH;
  int EangleH;
  int BangleL;
  int EangleL;
  int qL;    // quality of detection
  int qH;    // quality of detection
} BDL_IrisParams;

// verify eye parameters are sensible
EXTERNC MIA_RESULT_CODE BDL_VerifyEyeParams(
  const SPupilInfo* psPI,         // IN:  pupil parameters
  const SIrisInfo* psII,          // IN:  iris parameters
  int xs,                   // IN:  image size
  int ys);                  //

// generate flash mask, Seoul version
EXTERNC MIA_RESULT_CODE BDL_Flash_DetectFlash(
  unsigned char* pImage,
  int nWidth,
  int nHeight,
  unsigned char* pFlashMask);

EXTERNC MIA_RESULT_CODE BDL_Flash_RemoveFlashEx(
  unsigned char* pImage,
  int nWidth,
  int nHeight,
  unsigned char* pFlashMask);

EXTERNC MIA_RESULT_CODE BDL_Flash_RemoveFlash(
  unsigned char* pImage,
  int nWidth,
  int nHeight);

// generate flash mask, Moscow version
EXTERNC MIA_RESULT_CODE BDL_Flash_DetectMask_v2(
  unsigned char* mask,        // OUT: flash mask
  const unsigned char* im,    // IN:  input image
  int xs,                     // IN:  image size
  int ys,                     //
  const char* nambeg);

EXTERNC MIA_RESULT_CODE BDL_Flash_Enumerate(
  SConObjPrimaryParams** ppCOPP,    // OUT: objects
  int* pnCOPP,                      // OUT: number of objects
  const uint8* im,                  // IN:  source image
  int xs,                           // IN:  source size
  int ys,
  int thr);                         // IN:  brightness threshold for flash

// Center_Search function purpose is to find center of image (eye)
// the function returns (asCI, SCenterInfo is defined in bdl_struct.h) detected coordinates (xc,yc)
// and provides some quality of center detection asCI.q, which has not big meaning for this function
// 'mode' flag points the mode of center searching algorithm (defined in bdl_struct.h)
// IRL_CENTER_FLAG_CALC_FAST works quicker
// IRL_CENTER_FLAG_CALC_ROBUST works slower but more precise
EXTERNC MIA_RESULT_CODE BDL_Center_Search(
  SCenterInfo* asCI,          // OUT: center coordinates
  const unsigned char* im,    // IN:  input image
  int xs,                     // IN:  image size
  int ys,                     //
  int mode,                   // IN:  flags, see IRL_CENTER_xxx
  void* pvBuf,                // IN:  temporary buffer
  int* pnLen,                 // IN/OUT: allocated/used
  const char* nambeg);        // IN:  set this to NULL

//Bazrad_Search inputs preliminary detected center of eye and detects a circle of iris
//The function returns (in SBazRadInfo) coordinates of supposed iris circle parameters (x, y, r)
//Bazrad_Search function tries to detect rough pupil circle
//for that purpose SBazRadInfo structure contains pupil circle fields 
EXTERNC MIA_RESULT_CODE BDL_Bazrad_Search(
  SBazRadInfo* asBRI,             // OUT: eye BazRad evaluation
  const unsigned char* im,        // IN:  source image
  const SCenterInfo* pCI,        // IN: center coordinates
  int xs,                         // IN:  image size
  int ys,                         //
  int nBazradPredef,              // IN: predefined value of iris radius, -1 = no predefinition
  int nBazRadSpread,              // IN: spread around previous argument.
                                  //   i.e. radii may actually be Predef+-Spread
  void* buf,                      // IN:  buffer
  int* pnSz,                      // IN/OUT: size allocated/used
  const char* pcNamBeg);          // IN:  for debugging, set it to NULL

//Bazrad_SearchExt inputs preliminary detected center of eye and detects a circle of iris
//The function returns pnRads hypothesis (in SBazRadInfo array) set of supposed coordinates 
//of iris circle parameters (x, y, r) with their internal estimation
//Bazrad_SearchMany function tries to detect rough pupil circle too
//for that purpose SBazRadInfo structure contains pupil circle fields 
EXTERNC MIA_RESULT_CODE BDL_Bazrad_SearchMany(
  SBazRadInfo* asBRI,             // OUT: eye BazRad evaluation
  int* pnRads,                    // IN/OUT: number of positions allocated/used
  const unsigned char* im,        // IN:  source image
  const SCenterInfo* pCI,         // IN: center coordinates
  int xs,                         // IN:  image size
  int ys,                         //
  void* buf,                      // IN:  buffer
  int* pnSz,                      // IN/OUT: size allocated/used
  const char* pcNamBeg);          // IN:  for debugging, set it to NULL

// get pupil parameters in polar coordinates
//Pupil_Search inputs preliminary detected eye center and radius of iris (base radius)
//BDL_PupilParams output structure contains all computed in the function pupil details:
//pupil circle, pupil ellipse, etc.
EXTERNC MIA_RESULT_CODE BDL_PupilCspath_Search(
  const unsigned char* pBufImg,     // IN: eye image
  int  SizeX,                       // IN: image size
  int  SizeY,                       //
  int  xc,                          // IN: coordinates of eye Center 
  int  yc,                          //
  SBazRadInfo*  pBazRad,            // IN: coefficient for scale
  BDL_PupilParams* pPupParam,
  void*            MemoryBuff,      // IN:  tmp memory buffer
  int*             MemoryBuffSize); // IN:  tmp memory buffer size from GetParamPupil_GetMemBuffSize

#define IRL_GRAPAR_EDGEIMAGE    0x00000001  // edge image is passed
#define IRL_GRAPAR_GRADIMAGE    0x00000002  // gradient image is ready
#define IRL_GRAPAR_CALCRADIUS   0x00000004  // perform radius calculation

// search pupil
EXTERNC MIA_RESULT_CODE BDL_PupilGrapar_Search(
  SPupilInfo* sPI,    // OUT: pupil parameters
  const uint8* im,    // IN:  eye image
  const int16* pgx,   // IN:  gradient (if calculated)
  const int16* pgy,   // 
  int xs,             // IN:  image size
  int ys,             //
  int nFlags,         // IN:  processing parameters
  void* pvBuf,        // IN:  tmp memory buffer
  int*  pnLen,        // IN/OUT: bytes allocated/used
  const char* nambeg);// for debug, set it to NULL

// get pupil parameters in polar coordinates
//Pupil_Search inputs preliminary detected eye center and radius of iris (base radius)
//BDL_PupilParams output structure contains all computed in the function pupil details:
//pupil circle, pupil ellipse, etc.
EXTERNC MIA_RESULT_CODE BDL_PupilTriang_Search(
  const unsigned char* pBufImg,     // IN: eye image
  int  SizeX,                       // IN: image size
  int  SizeY,                       //
  int  xc,                          // IN: coordinates of eye Center 
  int  yc,                          //
  SBazRadInfo*  pBazRad,            // IN: coefficient for scale
  BDL_PupilParams* pPupParam,
  void*            MemoryBuff,      // IN:  tmp memory buffer
  int*             MemoryBuffSize); // IN:  tmp memory buffer size from GetParamPupil_GetMemBuffSize

//Iris_Search inputs preliminary detected eye pupil parameters and radius of iris (base radius)
//BDL_IrisParams output structure contains all computed in the function iris details:
//iris circle, iris ellipse, etc.
EXTERNC MIA_RESULT_CODE BDL_Iris_Search(
  const unsigned char *pBufImg, // IN: eye image
  int      SizeX,               // IN: image size
  int      SizeY,               //
  const    SPupilInfo* pPupilParam,
  const    SBazRadInfo* pBR,
  BDL_IrisParams*      pIrisParam,
  void*    MemoryBuff,      // IN:  tmp memory buffer
  int*     MemoryBuffSize); // IN:  tmp memory buffer size

EXTERNC MIA_RESULT_CODE IVIR_PupilFromIris_v1(
  SPupilInfo* pPupil,       // OUT: pupil data
  const SIrisInfo* pIris,   // IN:  iris data
  const unsigned char* im,  // IN:  source image
  int xs,                   // IN:  image size
  int ys,                   //
  void* buf,                // IN:  buffer
  int* buflen,              // IN/OUT: allocated/used buffer size
  const char* nambeg);      // IN:  used for debugging, set it to NULL

// get pupil parameters in polar coordinates
//Pupil_Search inputs preliminary detected eye center and radius of iris (base radius)
//BDL_PupilParams output structure contains all computed in the function pupil details:
//pupil circle, pupil ellipse, etc.
EXTERNC MIA_RESULT_CODE BDL_RefinePupil_Search(
  const unsigned char* pBufImg,     // IN: eye image
  int  SizeX,                       // IN: image size
  int  SizeY,                       //
  int  xc,                          // IN: coordinates of eye Center 
  int  yc,                          //
  BDL_PupilParams* pPupParam,
  void*            MemoryBuff,      // IN:  tmp memory buffer
  int*             MemoryBuffSize)  ;// IN:  tmp memory buffer size from GetParamPupil_GetMemBuffSize

EXTERNC MIA_RESULT_CODE BDL_RefinePupil_EllBlob_v1(
  SEllipseExt* pEll,    // OUT: ellipse with qualities
  uint8* im,      // IN:  source image
  int xs,               // IN:  image size
  int ys,
  const char* nambeg);  // IN:  for debuf, set it to NULL

EXTERNC MIA_RESULT_CODE BDL_RefinePupil_EllBlob_v2(
  SEllipseExt* pEll,    // OUT: ellipse with qualities
  uint8* im,      // IN:  source image
  int xs,               // IN:  image size
  int ys,
  int bCutFlash,        // IN:  do flash removal
  const char* nambeg);  // IN:  for debuf, set it to NULL

EXTERNC MIA_RESULT_CODE BDL_RefinePupil_EllBlob_v3(
  SEllipseExt* pEll,    // OUT: ellipse with qualities
  uint8* im,            // IN:  source image
  int xs,               // IN:  image size
  int ys,
  int bCutFlash,        // IN:  do flash removal
  const char* nambeg);  // IN:  for debuf, set it to NULL

//Occlusion_Search inputs preliminary detected eye iris parameters
//the function purpose is to make coherent region of iris suitable for recognition
//region is pointed by angles from iris center
EXTERNC MIA_RESULT_CODE BDL_RefineIris_Search(
  const unsigned char *bufin,  //IN: -input array
  int SizeX,            //IN: -size for X & Y for bufin
  int SizeY,            //IN:
  int *XIris_,          //IN:- parameters of Iris
  int *YIris_,          //IN:
  int *RIris_,          //IN:
  short *Buflids,       //IN: - array for lids
  int *BangleH,         //
  int *EangleH,
  int *BangleL,
  int *EangleL,
  int SizeBuf);         //IN: -size for Buflids and AR_POINT_L  

// iris border tracking by gradient
EXTERNC MIA_RESULT_CODE BDL_RefineIris_Occlusion_v5(
  unsigned char* mask,
  const unsigned char* im,    // IN:  input image
  int xs,                     // IN:  image size
  int ys,                     //
  const SPupilInfo* psPI,
  const SIrisInfo* psII,
  const char* nambeg);

// iris border tracking by gradient. same as _v5 but main part is scaled down by 2.
EXTERNC MIA_RESULT_CODE BDL_RefineIris_Occlusion_v6(
  unsigned char* mask,
  const unsigned char* im,    // IN:  input image
  int xs,                     // IN:  image size
  int ys,                     //
  const SPupilInfo* psPI,
  const SIrisInfo* psII,
  const char* nambeg);

MIA_RESULT_CODE BDL_RefineIris_CalculateOverlapping(
  int* pnQuality,
  const uint8* pucMask,
  int xs,
  int ys,
  const SPupilInfo* psPI,
  const SIrisInfo* psII);

// create image of the mask from polar representation of occluded/nonoccluded
EXTERNC MIA_RESULT_CODE BDL_RefineIris_CreatemaskFromBorderMarks(
  unsigned char* mask,    // OUT: mask image
  int xs,                 // IN:  image size
  int ys,                 //
  const int* goods,       // IN:  polar representation
  int len,                // IN:  number of points in polar
  const SPupilInfo* psPI, // IN:  pupil data
  const SIrisInfo* psII); // IN:  iris data

// exported here for debug / temporary
EXTERNC void BDL_RefineIris_EliminateOcclusionGaps(
  int* dst,
  int len);
// flags
#define BDL_LIVENESS_DITHER_FOURIER_VER1    0x00000101
#define BDL_LIVENESS_DITHER_FOURIER_VER2    0x00000102
#define BDL_LIVENESS_DITHER_GRADHIST_VER1   0x00000201
#define BDL_LIVENESS_DITHER_GRADHIST_VER2   0x00000202
#define BDL_LIVENESS_DITHER_GRADHIST_VER3   0x00000203
#define BDL_LIVENESS_DITHER_MORPH_VER1      0x00000301
#define BDL_LIVENESS_DITHER_MORPH_VER2      0x00000302
#define BDL_LIVENESS_DITHER_MORPH_VER3      0x00000303
#define BDL_LIVENESS_DITHER_DEFAULT         BDL_LIVENESS_DITHER_GRADHIST_VER3
#define BDL_LIVENESS_DITHER_FAST            BDL_LIVENESS_DITHER_GRADHIST_VER3
#define BDL_LIVENESS_DITHER_PRECISE         BDL_LIVENESS_DITHER_MORPH_VER3

//============ dither detection
EXTERNC MIA_RESULT_CODE BDL_Liveness_DitherDetect(
  double *pThr,         // OUT: fake probability
  const uint8* imgBuf,  // IN:  image buffer
  int xs,               // IN:  image width
  int ys,               // IN:  image height
  int nFlags,           // IN:  flags for method
  const char* nambeg);  // IN:  for debug, set to NULL

//============ pattern detection
DECL_HANDLE(HPatternAnalyzer);

typedef struct
{
  int nFlash1;                  // IN:  presence of flashes
  int nFlash2;
  int nFlash3;
  int nFlash4;
} SFlashInfo;

// create analyzer structure
EXTERNC MIA_RESULT_CODE BDL_Liveness_PatternCreateAnalyzer(
  HPatternAnalyzer* phAnalyzer);  // OUT: handler to pattern analyzer

// process one frame: detect flashes in image and take flash information from HW
EXTERNC MIA_RESULT_CODE BDL_Liveness_PatternProcessFrame(
  HPatternAnalyzer hAnalyzer,   // IN:  handler to pattern analyzer
  const unsigned char* im,      // IN:  image
  int xs,                       // IN:  image width
  int ys,                       // IN:  image height
  const SPupilInfo* pu,         // IN:  pupil data,
  const SIrisInfo* ir,          // IN:  iris data
  const SFlashOnInfo* fl,      // IN:  flash info (about turned on flashes)
  const char* nambeg);         // IN:  for debug, set it to NULL

// analyze accumulated sequence of frames
EXTERNC MIA_RESULT_CODE BDL_Liveness_PatternAnalyze(
  double* pdLiveness,           // OUT: liveness in range [0;1]
  HPatternAnalyzer hAnalyzer);  // IN:  handler to pattern analyzer

// free all resources of analyzer
EXTERNC MIA_RESULT_CODE BDL_Liveness_PatternFreeAnalyzer(
  HPatternAnalyzer* phAnalyzer);

//============ combined detection
DECL_HANDLE(HCombinedLivenessAnalyzer);

#define BDL_LIVENESS_COMBI_USEDITHER    0x00000001
#define BDL_LIVENESS_COMBI_USEPATTERN   0x00000002
#define BDL_LIVENESS_COMBI_DEFAULT      0x00000003

// create analyzer structure
EXTERNC MIA_RESULT_CODE BDL_Liveness_CombinedCreateAnalyzer(
  HCombinedLivenessAnalyzer* phAnalyzer,  // OUT: handler to analyzer
  int nFlags);                            // IN:  flags

// process one frame: detect flashes in image and take flash information from HW
EXTERNC MIA_RESULT_CODE BDL_Liveness_CombinedProcessFrame(
  HCombinedLivenessAnalyzer hAnalyzer,    // IN:  handler to pattern analyzer
  const unsigned char* im,      // IN:  image
  int xs,                       // IN:  image width
  int ys,                       // IN:  image height
  const SPupilInfo* pu,         // IN:  pupil data,
  const SIrisInfo* ir,          // IN:  iris data
  const SFlashOnInfo* fl);      // IN:  flash info (about turned on flashes)

// analyze accumulated sequence of frames
EXTERNC MIA_RESULT_CODE BDL_Liveness_CombinedAnalyze(
  double* pdLiveness,           // OUT: liveness in range [0;1]
  HCombinedLivenessAnalyzer hAnalyzer);   // IN:  handler to pattern analyzer

// free all resources of analyzer
EXTERNC MIA_RESULT_CODE BDL_Liveness_CombinedFreeAnalyzer(
  HCombinedLivenessAnalyzer* phAnalyzer);

// in fact, this is an internal function for occlusion detection
// inserted here temporary for debugging purposes
EXTERNC MIA_RESULT_CODE BDL_Occlusion_GetBestFrame_wrapper(
  int *pnBestPos,
  const unsigned char* im,
  int xs,
  int ys,
  const SPupilInfo* psPI,
  const SIrisInfo* psII,
  const char* nambeg);

// in fact, this is an internal function for occlusion detection
// inserted here temporary for debugging purposes
EXTERNC void BDL_Occlusion_GetBestFrame(
  int* pnFrameCenter,
  const unsigned char* im,
  int xs,
  int top,
  int bot,
  int wid,
  void* pvMem);  // xs*sizeof(int) bytes

// get mask of occlusion
MIA_RESULT_CODE BDL_Occlusion_Mask_v1(
  unsigned char* mask,        // OUT: mask
  const unsigned char* im,    // IN:  input image
  int xs,                     // IN:  image size
  int ys,                     //
  const SPupilInfo* psPI,     // IN:  pupil
  const SIrisInfo* psII,      // IN:  iris
  const char* nambeg);        // IN:  for debug, set it to NULL

#endif //__bdl_h__
