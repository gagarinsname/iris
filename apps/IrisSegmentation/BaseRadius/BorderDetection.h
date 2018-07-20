/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  23 April 2002                                          */
/*      Purpose:  Iris border detection API (common part                 */
/*      Authors:                                                         */
/*        Evgeny Demin                                                   */
/*        Konstantin Gankin                                              */
/*        Ivan Matveev                                                   */
/*        Iurii Efimov                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#ifndef __BorderDetection_h__
#define __BorderDetection_h__

#include "GlobalParams.h"

// #include "errorcodes.h"
// #include "stddefs.h"
// #include "ipl.h"

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
	int x, y;
} sPnt;

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
    unsigned char* TotalMask;
    unsigned char* FlashMask;
    unsigned char* OcclusionMask;
    int nTotalMaskQuality;  // percent of coverage of iris with mask
    int nDetQuality;
    // debug data
    char* acFileNameBeginning;
} SBorderData_ALL;


/*
// verify eye parameters are sensible
EXTERNC RESULT_CODE BD_VerifyEyeParams(
    const SPupilInfo* psPI,         // IN:  pupil parameters
    const SIrisInfo* psII,          // IN:  iris parameters
    int xs,                         // IN:  image size
    int ys);                        //

EXTERNC RESULT_CODE BD_Center_Search(
    SCenterInfo* asCI,          // OUT: center coordinates
    const unsigned char* im,    // IN:  input image
    int xs,                     // IN:  image size
    int ys,                     //
    int mode,                   // IN:  flags, see IRL_CENTER_xxx
    void* pvBuf,                // IN:  temporary buffer
    int* pnLen,                 // IN/OUT: allocated/used
    const char* nambeg);        // IN:  set this to NULL
*/
//Bazrad_Search inputs preliminary detected center of eye and detects a circle of iris
//The function returns (in SBazRadInfo) coordinates of supposed iris circle parameters (x, y, r)
//Bazrad_Search function tries to detect rough pupil circle
//for that purpose SBazRadInfo structure contains pupil circle fields 
/*EXTERNC MIA_RESULT_CODE BDL_Bazrad_Search(
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
    */

//Bazrad_SearchExt inputs preliminary detected center of eye and detects a circle of iris
//The function returns pnRads hypothesis (in SBazRadInfo array) set of supposed coordinates 
//of iris circle parameters (x, y, r) with their internal estimation
//Bazrad_SearchMany function tries to detect rough pupil circle too
//for that purpose SBazRadInfo structure contains pupil circle fields 
/*EXTERNC RESULT_CODE BDL_Bazrad_SearchMany(
    SBazRadInfo* asBRI,             // OUT: eye BazRad evaluation
    int* pnRads,                    // IN/OUT: number of positions allocated/used
    const unsigned char* im,        // IN:  source image
    const SCenterInfo* pCI,         // IN: center coordinates
    int xs,                         // IN:  image size
    int ys,                         //
    void* buf,                      // IN:  buffer
    int* pnSz,                      // IN/OUT: size allocated/used
    const char* pcNamBeg);          // IN:  for debugging, set it to NULL
*/


/*
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
*/

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

#endif //__BorderDetection_h__