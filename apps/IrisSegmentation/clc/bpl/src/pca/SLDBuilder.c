/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:     8 December 2004                                     */
/*      Revision:    1.0.00                                              */
/*      Purpose:     Linear Discriminant Building                        */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#include <math.h>
#include <memory.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include "SLDBuilder.h"
#include "bpl.h"

MIA_RESULT_CODE BPL_LDB_Allocate(
  SLDBuilder** ppSLDB,
  int nSamples,   // number of samples  N
  int nDim,       // dimension of space M
  int nExp)       // necessary vectors  L
{
  SLDBuilder* pSLDB;
  int covdim;

    // check arguments
    if (ppSLDB==NULL)
      return ERR_GEN_NULLPOINTER;
    if ((nSamples==0)||(nDim==0)||(nExp==0))
      return ERR_GEN_INVALID_PARAMETER;
    if ((nExp>nSamples)||(nExp>nDim))
      return ERR_GEN_INVALID_PARAMETER;
    // decide the correlation method
    covdim = min(nDim,nSamples);
    // allocate
    if ((pSLDB = (SLDBuilder*)malloc(
          sizeof(SLDBuilder)+           // self
          nSamples*sizeof(int)+         // class indexes                N
          nSamples*nDim*sizeof(double)+ // samples in original space    N*M
//          nClass*nDim*sizeof(double)+   // mean sample vectors          K*M
//          nDim*nDim*sizeof(double)+     // inside-class covariation     M*M
//          nDim*nDim*sizeof(double)+     // between-class covariation    M*M
          nDim*nDim*sizeof(double)+     // eigenvectors of matrix       M*M
          nDim*sizeof(double)))==NULL)  // eigen values of matrix       M
      return ERR_GEN_NOMEMORY;
    // set pointers
    pSLDB->pnClassIndexes = (int*)(pSLDB+1);
    pSLDB->Samples = (double*)(pSLDB->pnClassIndexes+nSamples);
//    pSLDB->MeanSamples = pSLDB->Samples+nSamples*nDim;
//    pSLDB->InnerCov = pSLDB->MeanSamples+nClass*nDim;
//    pSLDB->OuterCov = pSLDB->InnerCov+nDim*nDim;
    pSLDB->EigenVec = pSLDB->Samples+nSamples*nDim;//pSLDB->OuterCov+nDim*nDim;
    pSLDB->EigenVal = pSLDB->EigenVec+nDim*nDim;
    // set data
    memset(pSLDB->pnClassIndexes,-1,nSamples*sizeof(pSLDB->pnClassIndexes[0]));
    pSLDB->nSampleCount = nSamples;
    pSLDB->nSampleDimension = nDim;
    pSLDB->nNecessaryVectors = nExp;
    // output
    *ppSLDB = pSLDB;
    return ERR_OK;
}

MIA_RESULT_CODE BPL_LDB_PutSourceImage(
  SLDBuilder* pSLDB,
  const double* im,   // image itself
  int num,            // position in source sequence
  int cls)            // class index
{
    // check arguments
    if ((pSLDB==NULL)||(im==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((num<0)||(num>=pSLDB->nSampleCount)||(cls<0))
      return ERR_GEN_INVALID_PARAMETER;
    // set data
    memcpy(pSLDB->Samples+num*pSLDB->nSampleDimension,im,pSLDB->nSampleDimension*sizeof(pSLDB->Samples[0]));
    pSLDB->pnClassIndexes[num] = cls;
    return ERR_OK;
}

//staticFunc int __cdecl ownLDB_SortAscending_CompareInt(
staticFunc int ownLDB_SortAscending_CompareInt(
  const int* el1,
  const int* el2)
{
    return *el1 - *el2;
}

// sort Values with corresponding Elements in value decreasing order
staticFunc MIA_RESULT_CODE ownSLDB_SortClassAscending(
  int* pValues,     // IN: values which give the sort order
  void* pElements,  // IN: elements
  int nNumber,      // IN: number of elements
  int nSize,        // IN: size of each element in bytes
  void* extbuftmp,  // IN: temporary buffer of size
  int nTMP)         // IN: size of allocated memory
{
  typedef struct
  {
    int val;
    int idx;
  } IDX_AND_VAL;
  //typedef int (__cdecl* cmp_fn_void)(const void *elem1, const void *elem2);
  typedef int(* cmp_fn_void)(const void *elem1, const void *elem2);

  IDX_AND_VAL* iav;
  char* elstore;
  int i,empty,tmpint;

    // check space
    if (nTMP<(int)(sizeof(IDX_AND_VAL)*nNumber+nSize))
      return ERR_GEN_NOMEMORY;
    // allocate memory for index<->value array
    iav = (IDX_AND_VAL*)extbuftmp;
    elstore = (char*)(iav+nNumber);
    // fill the indexing array
    for (i=0;i<nNumber;i++)
    {
      iav[i].val = pValues[i];
      iav[i].idx = i;
    }
    // sort the indexing array
    qsort(iav,nNumber,sizeof(IDX_AND_VAL),(cmp_fn_void)(&ownLDB_SortAscending_CompareInt));
    // output sorted sequence
    for (i=0;i<nNumber;i++)
      pValues[i] = iav[i].val;
    // rearrange element by the results of sorting
    // search for mismatch
    for (i=0;i<nNumber;i++)
      if (iav[i].idx!=i)
      { // mismatch found!
        // copy to a temporary safe
        memcpy(elstore,((char*)pElements)+i*nSize,nSize);
        // i-th element is emptied
        empty = i;
        // cycle
        while (iav[empty].idx!=i)
        {
          memcpy( ((char*)pElements)+empty*nSize,
                  ((char*)pElements)+iav[empty].idx*nSize,
                  nSize);
          tmpint = empty;
          empty = iav[empty].idx;
          iav[tmpint].idx = tmpint;
        }
        memcpy(((char*)pElements)+empty*nSize,elstore,nSize);
        iav[empty].idx = empty;
      }
    return ERR_OK;
}

// calculate total per-element covariation
staticFunc void ownSLDB_CovariateTotal(
  double* totcov,
  const double* Samples,
  int nSampleDimension,
  int nSampleCount)
{
  int i,j;
  const double *ptr,*finish_ptr;
  double sum_xy,sum_x,sum_y;

    for (i=0;i<nSampleDimension;i++)
      for (j=0;j<=i;j++)
      {
        sum_x = sum_y = sum_xy = 0.;
        ptr = Samples+j;
        finish_ptr = ptr+nSampleCount*nSampleDimension;
        for (;ptr != finish_ptr;ptr+=nSampleDimension)  // will do 'samplecount' loops
        {
          sum_x += ptr[0];
          sum_y += ptr[i-j];
          sum_xy += ptr[0]*ptr[i-j];
        }
        totcov[i*nSampleDimension+j] = totcov[j*nSampleDimension+i] =
          nSampleCount*sum_xy-sum_x*sum_y;
      }
}

// calculate the covariation of samples
staticFunc void ownSLDB_CovariateInner(
  double* inncov,    // OUT: size = nSampleDimension*nSampleDimension
  const double* Samples,
  int nSampleDimension,
  int nSampleCount,
  int nClassCount,
  const int *pnClassMemberCount)
{
  int i,j,k;
  const double *ptr,*finish_ptr;
  double sum_xy,sum_x,sum_y,sum_var;

    for (i=0;i<nSampleDimension;i++)
      for (j=0;j<=i;j++)
      {
        sum_var = 0;
        ptr = Samples+j;
        for (k=0;k<nClassCount;k++)
        {
          sum_x = sum_y = sum_xy = 0.;
          finish_ptr = ptr+pnClassMemberCount[k]*nSampleDimension;
          for (;ptr != finish_ptr;ptr+=nSampleDimension)  // will do 'ClassMemberCount[k]' loops
          {
            sum_x += ptr[0];
            sum_y += ptr[i-j];
            sum_xy += ptr[0]*ptr[i-j];
          }
          sum_var += pnClassMemberCount[k]*sum_xy-sum_x*sum_y;
        }
        inncov[i*nSampleDimension+j] = inncov[j*nSampleDimension+i] = sum_var;
      }
}

staticFunc void ownBPL_fnFakeLDCalculateCallback(SPCCalculateCallbackData* ptr) {};

// size of temporary buffer ismax of:
//    nSampleDimension*sizeof(int)
//    3*sizeof(double)*nSampleCount
MIA_RESULT_CODE BPL_LDB_Calculate(
  SLDBuilder* pSLDB,            // IN:  object
  pfnPCCalculateCallback pCB,   // IN:  callback function
  void* pvToken)                // IN:  caller's token
{
  MIA_RESULT_CODE res = ERR_OK;
  SPCCalculateCallbackData spccb; // passed to caller
  int i,nClassesCount,nCurrentClass,nBiggestClass,nSmallestClass,j;
  double ErrMaxAbsDevDiag;   // maximum absolute deviation of diagonal elements from 1
  double ErrMaxAbsDevNdiag;  // maximum absolute deviation of off-diagonal elements from 0
  double ErrAvgAbsDevDiag;   // average absolute deviation of diagonal elements from 1
  double ErrAvgAbsDevNdiag;  // average absolute deviation of off-diagonal elements from 0
  double ErrMaxDevEigen;     // maximum absolute deviation of diagonal elements from 1
  int64 tix;
  // work arrays
  double *InnerCov,*InnerCov_copy;  // inside-class covariation             M*M
  double *OuterCov,*OuterCov_copy;  // between-class covariation            M*M
  int *pnClassMemberCount;
  double tmpdbl;
  void* extbuftmp,*bup = NULL;  // temporary buffer
  int nTMP;                     // size of temporary buffer

    // check arguments
    if (pSLDB==NULL)
      return ERR_GEN_NULLPOINTER;
    // allocate work arrays
    if ((bup = malloc(
      4*pSLDB->nSampleDimension*pSLDB->nSampleDimension*sizeof(double)+
      (nTMP = 2*pSLDB->nSampleDimension*sizeof(double)+
              2*pSLDB->nSampleCount*sizeof(int))))==NULL)
      return ERR_GEN_NOMEMORY;
    InnerCov = (double*)bup;         // inside-class covariation             M*M
    OuterCov = InnerCov+pSLDB->nSampleDimension*pSLDB->nSampleDimension;         // between-class covariation            M*M
    InnerCov_copy = OuterCov+pSLDB->nSampleDimension*pSLDB->nSampleDimension;
    OuterCov_copy = InnerCov_copy+pSLDB->nSampleDimension*pSLDB->nSampleDimension;
    extbuftmp = OuterCov_copy+pSLDB->nSampleDimension*pSLDB->nSampleDimension;
    // make fake if necessary
    if (pCB==NULL)
      pCB = ownBPL_fnFakeLDCalculateCallback;
    // pass callback token to callback structure
    spccb.pvToken = pvToken;
    // calculate mean vector and normalise image set to a zero mean
    strcpy(spccb.repstr,"LDCalc: Sorting vectors by class numbers");
    (*pCB)(&spccb);
    tix = clock();
    // sort Values with corresponding Elements in value decreasing order
    if ((res = ownSLDB_SortClassAscending(
      pSLDB->pnClassIndexes,     // IN: values which give the sort order
      pSLDB->Samples,  // IN: elements
      pSLDB->nSampleCount,      // IN: number of elements
      pSLDB->nSampleDimension*sizeof(pSLDB->Samples[0]), // IN: size in bytes
      extbuftmp,  // IN: temporary buffe: sz = sampledim*sz(double)
      nTMP))!=ERR_OK)   // IN: size of allocated memory
      goto endldbcalc;
    tix = clock()-tix;
    sprintf(spccb.repstr,"  Done in %d clox",(int)tix);
    (*pCB)(&spccb);
    // check first class is 0 (will be -1, if some samples were skipped)
    if (pSLDB->pnClassIndexes[0]!=0)
    {
      res = ERR_GEN_NOT_INITIALISED;
      goto endldbcalc;
    }
    // count of classes should be a number of last class+1
    nClassesCount = pSLDB->pnClassIndexes[pSLDB->nSampleCount-1]+1;
    // should be at least two classes
    if (nClassesCount<2)
    {
      res = MIA_LDA_TOO_FEW_CLASSES;
      goto endldbcalc;
    }
    // allocate memory for numbers of class
    if (nTMP<(i = nClassesCount*sizeof(int)))  // classcounts
    {
      res = ERR_GEN_NOMEMORY;
      goto endldbcalc;
    }
    pnClassMemberCount = (int*)extbuftmp;
    extbuftmp = (void*)(((ptr_t)(extbuftmp))+i);
    nTMP -= i;
    memset(pnClassMemberCount,0,nClassesCount*sizeof(pnClassMemberCount[0]));
    // calculate number of samples in each class, check contiguity of class numbers
    nCurrentClass = 0;
    nBiggestClass = nSmallestClass = 0;
    for (i=0;i<pSLDB->nSampleCount;i++)
      if (pSLDB->pnClassIndexes[i]==nCurrentClass)
        pnClassMemberCount[nCurrentClass]++;
      else
      { // class is changing, first check the class we have just finished
        // check if class is adequately populated
        if (pnClassMemberCount[nCurrentClass]<3)
        {
          res = MIA_LDA_TOO_SMALL_CLASS;
          goto endldbcalc;
        }
        // update smallest and biggest class numbers
        if (pnClassMemberCount[nBiggestClass]<pnClassMemberCount[nCurrentClass])
          nBiggestClass = nCurrentClass;
        else
          if (pnClassMemberCount[nSmallestClass]>pnClassMemberCount[nCurrentClass])
            nSmallestClass = nCurrentClass;
        // check contiguity (i.e. if class indexes go with unit increment)
        if (pSLDB->pnClassIndexes[i]-nCurrentClass!=1)
        {
          res = MIA_LDA_UNCONTIGUIOS_INDEXES;
          goto endldbcalc;
        }
        // change the class
        pnClassMemberCount[++nCurrentClass]++;
      }
    // last one was left unchecked. verify it has same class as previous
    if (pSLDB->pnClassIndexes[pSLDB->nSampleCount-1]!=pSLDB->pnClassIndexes[pSLDB->nSampleCount-2])
    {
      res = MIA_LDA_TOO_SMALL_CLASS;
      goto endldbcalc;
    }
    // report
    sprintf(spccb.repstr,"LDCalc: Biggest class %d has %d members, smallest class %d has %d members",
      nBiggestClass,pnClassMemberCount[nBiggestClass],nSmallestClass,pnClassMemberCount[nSmallestClass]);
    (*pCB)(&spccb);
    // calculte inner-class covariation
    strcpy(spccb.repstr,"LDCalc: Calculating inner-class covariation");
    (*pCB)(&spccb);
    tix = clock();
    ownSLDB_CovariateInner(
      InnerCov,
      pSLDB->Samples,
      pSLDB->nSampleDimension,
      pSLDB->nSampleCount,
      nClassesCount,
      pnClassMemberCount);
    sprintf(spccb.repstr,"  Done in %d clox",(int)((clock()-tix)));
    (*pCB)(&spccb);
    // calculate total covariation
    strcpy(spccb.repstr,"LDCalc: Calculating total covariation");
    (*pCB)(&spccb);
    tix = clock();
    ownSLDB_CovariateTotal(
      OuterCov,
      pSLDB->Samples,
      pSLDB->nSampleDimension,
      pSLDB->nSampleCount);
    sprintf(spccb.repstr,"  Done in %d clox",(int)(clock()-tix));
    (*pCB)(&spccb);
    // subtract inner from total obtaining outer
    for (i=0;i<pSLDB->nSampleDimension*pSLDB->nSampleDimension;i++)
      OuterCov[i] -= InnerCov[i];
    // store cov matrixes
    memcpy(InnerCov_copy,InnerCov,pSLDB->nSampleDimension*pSLDB->nSampleDimension*sizeof(InnerCov[0]));
    memcpy(OuterCov_copy,OuterCov,pSLDB->nSampleDimension*pSLDB->nSampleDimension*sizeof(OuterCov[0]));
    // solve eigen problem
    strcpy(spccb.repstr,"LDCalc: Solving generalized eigen problem");
    (*pCB)(&spccb);
    tix = clock();
    if ((res = BPL_Matrix_EigenGeneralized_Symm(
      OuterCov,        // IN:  left-size matrix (A)
      InnerCov,        // IN:  right-side matrix (B)
      pSLDB->EigenVec,        // OUT: eigen vectors
      pSLDB->EigenVal,       // OUT: eigen values
      pSLDB->nSampleDimension,            // IN:  dimensionality
      extbuftmp,  // IN:  temporary buffer, sz = 2*sampledim*sz(double)
      nTMP))!=ERR_OK)         // IN:  size of temporary buffer
      goto endldbcalc;
    sprintf(spccb.repstr,"  Done in %d clox",(int)(clock()-tix));
    (*pCB)(&spccb);
    // vectors are in eigen-increase order. reverse the order
    for (i=0;i<pSLDB->nSampleDimension/2;i++)
    {
      // exchange vectors
      memcpy(extbuftmp,
             pSLDB->EigenVec+pSLDB->nSampleDimension*i,
             pSLDB->nSampleDimension*sizeof(pSLDB->EigenVec[0]));
      memcpy(pSLDB->EigenVec+pSLDB->nSampleDimension*i,
             pSLDB->EigenVec+pSLDB->nSampleDimension*(pSLDB->nSampleDimension-i-1),
             pSLDB->nSampleDimension*sizeof(pSLDB->EigenVec[0]));
      memcpy(pSLDB->EigenVec+pSLDB->nSampleDimension*(pSLDB->nSampleDimension-i-1),
             extbuftmp,
             pSLDB->nSampleDimension*sizeof(pSLDB->EigenVec[0]));
      // exchange values
      ErrMaxAbsDevDiag = pSLDB->EigenVal[i];
      pSLDB->EigenVal[i] = pSLDB->EigenVal[pSLDB->nSampleDimension-i-1];
      pSLDB->EigenVal[pSLDB->nSampleDimension-i-1] = ErrMaxAbsDevDiag;
    }
    // calc eigen errors
    strcpy(spccb.repstr,"LDCalc: checking eigen errors");
    (*pCB)(&spccb);
    // test AQAt=E
    BPL_Matrix_EstimateErrorNonUnitarGeneral(
      pSLDB->EigenVec,InnerCov_copy,pSLDB->nSampleDimension,
      &ErrMaxAbsDevDiag,&ErrMaxAbsDevNdiag,
      &ErrAvgAbsDevDiag,&ErrAvgAbsDevNdiag);
    sprintf(spccb.repstr,"  Eigen matrix non-unitarity:\n\tMaxDiag = %le,\n\tMaxNdiag = %le,\n\tAvgDiag = %le,\n\tAvgNdiag = %le",
      ErrMaxAbsDevDiag,ErrMaxAbsDevNdiag,ErrAvgAbsDevDiag,ErrAvgAbsDevNdiag);
    (*pCB)(&spccb);
    // normalize to unit
    for (i=0;i<pSLDB->nSampleDimension;i++)
    {
      tmpdbl = 0.;
      for (j=0;j<pSLDB->nSampleDimension;j++)
        tmpdbl += pSLDB->EigenVec[i*pSLDB->nSampleDimension+j]*pSLDB->EigenVec[i*pSLDB->nSampleDimension+j];
      if (tmpdbl==0.)
      {
        sprintf(spccb.repstr,"LDCalc: vector %d is zero",i);
        (*pCB)(&spccb);
        res = MIA_LIN_MATRIX_DEGENERATE;
        goto endldbcalc;
      }
      tmpdbl = sqrt(tmpdbl);
      for (j=0;j<pSLDB->nSampleDimension;j++)
        pSLDB->EigenVec[i*pSLDB->nSampleDimension+j] /= tmpdbl;
    }
    // test Ax=LBx
    BPL_Matrix_EstimateErrorNonEigenGeneral(
      OuterCov_copy,    // IN:  source matrix at left
      InnerCov_copy,    // IN:  source matrix at right
      pSLDB->EigenVec,  // IN:  eigen vectors
      pSLDB->EigenVal,  // IN:  eigen values
      pSLDB->nSampleDimension,              // IN:  size of source matrix
      pSLDB->nNecessaryVectors,              // IN:  number of eigen vectors to be tested
      &ErrMaxDevEigen);
    sprintf(spccb.repstr,"  Maximum of ||Ax-Lx||/L = %le",ErrMaxDevEigen);
    (*pCB)(&spccb);
    // ok
endldbcalc:
    if (bup)
      free(bup);
    return res;
}
