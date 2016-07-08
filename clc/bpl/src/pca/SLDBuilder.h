/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:     8 December 2004                                     */
/*      Revision:    1.0.00                                              */
/*      Purpose:     Linear Discriminant Building                        */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#ifndef __SLDBUILDER_H__
#define __SLDBUILDER_H__

#include "SPCBuilder.h"

typedef struct
{
// dimensions
  int nSampleCount;           // number of samples                  =N
  int nSampleDimension;       // dimension of samples               =M
          // a.k.a. number of parameters
          // if done over PCA then in real usage it should be       M<N
  int nNecessaryVectors;    // number of necessary vectors          =L
          // for majority of purposes we need                       L<<N
          // by method definition                                   L<=M
// source data
  int* pnClassIndexes;      // class index to which 
                            // certain sample belongs to            N
                            // also used as sample load stamps
  double* Samples;          // samples in original space            N*M
// resulting data
  double* EigenVec;         // eigenvectors of matrix               M*M
  double* EigenVal;         // eigen values of matrix               M
} SLDBuilder;

EXTERNC MIA_RESULT_CODE BPL_LDB_Allocate(
  SLDBuilder** ppSLDB,
  int nSamples,   // number of samples  N
  int nDim,       // dimension of space M
  int nExp);      // necessary vectors  L

#define SLDBuilder_Free(pSLDB) \
{ \
    if (pSLDB) \
      free(pSLDB); \
}

EXTERNC MIA_RESULT_CODE BPL_LDB_PutSourceImage(
  SLDBuilder* pSLDB,
  const double* im,   // image itself
  int num,            // position in source sequence
  int cls);           // class index

EXTERNC MIA_RESULT_CODE BPL_LDB_Calculate(
  SLDBuilder* pSLDB,          // IN:  object
  pfnPCCalculateCallback pCB, // IN:  callback function
  void* pvUD);                // IN:  callback user data

#endif //__SLDBUILDER_H__
