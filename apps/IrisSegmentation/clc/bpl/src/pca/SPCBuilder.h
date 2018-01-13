/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:     21 March 2001                                       */
/*      Modified:    10 September 2001 - add PCA expansion errors        */
/*      Modified:    11 May 2004 - CPP->C, double->int (v3)              */
/*      Revision:    3.0.00                                              */
/*      Purpose:     Principal Components Building                       */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#ifndef __SPCBUILDER_H__
#define __SPCBUILDER_H__

#include "../ownbpl.h"

typedef struct
{
// dimensions
  int nSampleCount;           // number of samples            =N
                              // i.e. dimension of full eigen space
  int nSampleDimension;       // dimension of samples         =M
                              // i.e. source space dimension
    // almost everywhere in image processing                 M>>N
  int nNecessaryVectors;      // number of necessary vectors  =L
    // for majority of purposes we need                      L<<N
  int isPerElem;              // what type of correlation is done
// source statistical data
  char* Samples;              // samples in original space    N*M
  int* Samples_Present;       // sample load stamps
// resulting data
  unsigned char* MeanSample;  // mean sample vector           M
  double* EigenSpace;         // eigenspace                   L*M
} SPCBuilder;

// structure coming from PCB_Calculate in a callback to caller's code
typedef struct
{
  void* pvToken;      // token from caller passed thru BuildRules arg
  char repstr[256];
} SPCCalculateCallbackData;

typedef void (*pfnPCCalculateCallback) (SPCCalculateCallbackData*);

EXTERNC MIA_RESULT_CODE BPL_PCB_Allocate(
  SPCBuilder** ppSPCB,
  int nSamples,   // number of samples N
  int nDim,       // dimension of space M
  int nExp);      // necessary vectors L

#define SPCBuilder_Free(pSPCB) \
{ \
    if (pSPCB) \
      free(pSPCB); \
}

EXTERNC MIA_RESULT_CODE BPL_PCB_PutSourceImage(
  SPCBuilder* pSPCB,
  const unsigned char* im,  // image itself
  int num);       // position in source sequence

#define BPL_PCB_GetSourceImage(pSPCB,n) ((pSPCB)->Samples+(n)*(pSPCB)->nSampleDimension)

EXTERNC MIA_RESULT_CODE BPL_PCB_Calculate(
  SPCBuilder* pSPCB,          // IN:  object
  pfnPCCalculateCallback pCB, // IN:  callback function
  void* pvUD);                // IN:  callback user data
/*
#define PC_BLOCK_SIZE(IDim,PDim) \
    (4+2+2+(IDim)+(SHORTSIZE/8)*(IDim)*(PDim)+(PDim))

#define PC_BLOCK_SIZE_ALIGNED(IDim,PDim) \
    PC_BLOCK_SIZE((((IDim)+SUMCHUNK-1)/SUMCHUNK)*SUMCHUNK,\
                  (((PDim)+PC_NVEC_ALIGN-1)/PC_NVEC_ALIGN)*PC_NVEC_ALIGN)

#define BPL_PCB_SavePC_SIZELOC(pSPCB) \
  PC_BLOCK_SIZE_ALIGNED(pSPCB->nSampleDimension,pSPCB->nSampleCount)
*/
// size of buffer is: PDim
EXTERNC MIA_RESULT_CODE BPL_PCB_SavePC(
  const SPCBuilder* pSPCB,      // IN: PC space to be saved
  unsigned char* dst,           // IN: buffer accepting space
  int* dstsize,                 // IN: size of buffer, OUT: bytes actually placed in buffer
  int PDim);                    // IN: space dimension

EXTERNC MIA_RESULT_CODE BPL_PCB_Expand(
  const SPCBuilder* pSPCB,  // IN:  PCA structure
  double* coefs,            // OUT: expansion coefficients
  const unsigned char* im); // IN:  source image

#endif //__SPCBUILDER_H__
