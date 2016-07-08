/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:     30 March 2001                                       */
/*      Modified:    5 June 2004 - PCanalyser and PCBuilder splitted (v2)*/
/*      Modified:    17 November 2004 - revised normalise and covar fns  */
/*      Modified:    18 November 2004 - added per-element covariation    */
/*      Revision:    2.2.00                                              */
/*      Purpose:     internal PCA class functions for buiding PC space   */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/
#define __FILENUM__ 97 // __FILENUM__TAG97

#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <time.h>
#include "SPCBuilder.h"
#include "bpl.h"

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))

MIA_RESULT_CODE BPL_PCB_Allocate(
  SPCBuilder** ppSPCB,
  int nSamples,   // number of samples N
  int nDim,       // dimension of space M
  int nExp)       // necessary vectors L
{
  SPCBuilder* pSPCB;
  int covdim;

    // check arguments
    if (ppSPCB==NULL)
      return ERR_GEN_NULLPOINTER;
    if ((nSamples==0)||(nDim==0)||(nExp==0))
      return ERR_GEN_INVALID_PARAMETER;
    if ((nExp>nSamples)||(nExp>nDim))
      return ERR_GEN_INVALID_PARAMETER;
    // decide the correlation method
    covdim = min(nDim,nSamples);
    // allocate
    if ((pSPCB = (SPCBuilder*)malloc(
          sizeof(SPCBuilder)+               // self
          nSamples*nDim*sizeof(char)+       // samples
          nSamples*sizeof(int)+             // sample present
          nDim*sizeof(unsigned char)+       // mean sample
//          covdim*covdim*sizeof(double)+     // covariation matrix
  //        covdim*covdim*sizeof(double)+     // eigen vectors
    //      covdim*sizeof(double)+            // eigen values
          nExp*nDim*sizeof(double)          // eigenspace
         ))==NULL)
      return ERR_GEN_NOMEMORY;
    // set pointers
    pSPCB->Samples = (char*)(pSPCB+1);
    pSPCB->Samples_Present = (int*)(pSPCB->Samples+nSamples*nDim);
    pSPCB->MeanSample = (unsigned char*)(pSPCB->Samples_Present+nSamples);
//    pSPCB->CovMatrix = (double*)(pSPCB->MeanSample+nDim);
  //  pSPCB->EigenVec = pSPCB->CovMatrix+covdim*covdim;
    //pSPCB->EigenVal = pSPCB->EigenVec+covdim*covdim;
    pSPCB->EigenSpace = (double*)(pSPCB->MeanSample+nDim);//pSPCB->EigenVal+covdim;
    // set data
    memset(pSPCB->Samples_Present,0,nSamples*sizeof(pSPCB->Samples_Present[0]));
    pSPCB->nSampleCount = nSamples;
    pSPCB->nSampleDimension = nDim;
    pSPCB->nNecessaryVectors = nExp;
    pSPCB->isPerElem = (nSamples>=nDim)?1:0;;
    // output
    *ppSPCB = pSPCB;
    return ERR_OK;
}

MIA_RESULT_CODE BPL_PCB_PutSourceImage(
  SPCBuilder* pSPCB,
  const unsigned char* im,  // image itself
  int num)        // position in source sequence
{
    // check arguments
    if ((pSPCB==NULL)||(im==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((num<0)||(num>=pSPCB->nSampleCount))
      return ERR_GEN_INVALID_PARAMETER;
    // set data
    memcpy(pSPCB->Samples+num*pSPCB->nSampleDimension,im,pSPCB->nSampleDimension);
    pSPCB->Samples_Present[num]++;
    return ERR_OK;
}

staticFunc MIA_RESULT_CODE ownSPCB_NormaliseOnMean(
  unsigned char* MeanSample,
  const char* Samples,
  int nSampleDimension,
  int nSampleCount,
  void* extbuftmp,      // temporary buffer of size nSampleDimension*sizeof(int)
  int nTMP)
{
  unsigned int* Mean;
  int i,j,val;
  unsigned char* ptr;

    // allocate
    if (nTMP<nSampleDimension*(int)sizeof(Mean[0]))
      return ERR_GEN_NOMEMORY;
    Mean = (unsigned int*)extbuftmp;
    memset(Mean,0,nSampleDimension*sizeof(Mean[0]));
    // sum
    for (j=0;j<nSampleCount;j++)
    {
      ptr = (unsigned char*)(Samples+j*nSampleDimension);
      for (i=0;i<nSampleDimension;i++)
        Mean[i] += ptr[i];
    }
    // divide
    for (i=0;i<nSampleDimension;i++)
      MeanSample[i] = (unsigned char)((Mean[i]+nSampleCount/2)/nSampleCount);
    // normalise source images
    for (j=0;j<nSampleCount;j++)
    {
      ptr = (unsigned char*)(Samples+j*nSampleDimension);
      for (i=0;i<nSampleDimension;i++)
      {
        val = (int)((unsigned int)(ptr[i]))-(int)((unsigned int)(MeanSample[i]));
        if (val>=-127)  // I HATE this integer 0x80XX asymmetry! better make this NaN
        {
          if (val<127)
            ((char*)ptr)[i] = (char)val;
          else
            ((char*)ptr)[i] = 127;
        }
        else
          ((char*)ptr)[i] = -127;
      }
    }
    return ERR_OK;
}

// calculate the covariation of samples
staticFunc void ownSPCB_CovariatePerSamples(
  double* SampleCov,
  const char* Samples,
  int nSampleDimension,
  int nSampleCount)
{
  int64 counter64; // here the sum is limited by 2^36 so we use int64
  int i,j,offset;
  const char *ptr,*finish_ptr;

    // perform the calculation of uncentered moments
    for ( i=0;i<nSampleCount;i++)
      for (j=0;j<=i;j++)
      {
        counter64 = 0;
        ptr = Samples+j*nSampleDimension;
        finish_ptr = ptr+nSampleDimension;
        offset = nSampleDimension*(i-j);
        for (;ptr != finish_ptr;ptr++)
          counter64 += ((short)ptr[0])*((short)ptr[offset]);
        SampleCov[i*nSampleCount+j] = SampleCov[j*nSampleCount+i] = ((double)counter64)/nSampleDimension;
      }
}

// calculate the covariation of samples
staticFunc void ownSPCB_CovariatePerElements(
  double* ElemCov,    // OUT: size = nSampleDimension*nSampleDimension
  const char* Samples,
  int nSampleDimension,
  int nSampleCount)
{
  int64 counter64; // here the sum is limited by 2^36 so we use int64
  int i,j,offset;
  const char *ptr,*finish_ptr;

    // perform the calculation of uncentered moments
    for ( i=0;i<nSampleDimension;i++)
      for (j=0;j<=i;j++)
      {
        counter64 = 0;
        ptr = Samples+j;
        finish_ptr = ptr+nSampleCount*nSampleDimension;
        offset = i-j;
        for (;ptr != finish_ptr;ptr+=nSampleDimension)
          counter64 += ((short)ptr[0])*((short)ptr[offset]);
        ElemCov[i*nSampleDimension+j] = ElemCov[j*nSampleDimension+i] = ((double)counter64)/nSampleCount;
      }
}

// multiplying double matrix by char matrix
// used for calculation of F=B*V*L=Vt*Bt*L
// where V - eigen vectors, B - samples, L - inverse sqrt of eigen values
staticFunc void ownSPCB_MatrixSpecProductABL(
  double* C,          // destination: double matrix IxK
  double* A,    // source1: double matrix     IxJ
  const char* B,      // source2: char matrix       JxK
  const double* v,    // source3: values            I
  int I,              // size I   (number of necesseary eigen vectors)
  int J,              // size J   (number of samples, J>=I)
  int K)              // size K   (dimensionality of samples, K>>J)
{
  int i,j,k;
  double sum,val;

    // multiple small matrices V and L
    // this is faster than multiplying by L in major cycle
    for (i=0;i<I;i++)
    {
      val = 1./sqrt(v[i]*K);
      for (j=0;j<J;j++)
        A[i*J+j] *= val;
    }
    // major cycle
    for (i=0;i<I;i++)
    {
      for (k=0;k<K;k++)
      {
        sum = 0.;
        for (j=0;j<J;j++)
          sum += A[i*J+j]*B[j*K+k];
        C[i*K+k] = sum;
      }
    }
}

staticFunc void ownBPL_fnFakePCCalculateCallback(SPCCalculateCallbackData* ptr) {};

typedef void (*pfnCovariation)(double*,const char*,int,int);

MIA_RESULT_CODE BPL_PCB_Calculate(
  SPCBuilder* pSPCB,        // IN:  object
  pfnPCCalculateCallback pCB,       // IN:  callback function
  void* pvToken)            // IN:  caller's token
{
  MIA_RESULT_CODE res = ERR_OK;
  int i,startclear,eigendim;
  double ErrMaxAbsDevDiag;   // maximum absolute deviation of diagonal elements from 1
  double ErrMaxAbsDevNdiag;  // maximum absolute deviation of off-diagonal elements from 0
  double ErrAvgAbsDevDiag;   // average absolute deviation of diagonal elements from 1
  double ErrAvgAbsDevNdiag;  // average absolute deviation of off-diagonal elements from 0
  double ErrMaxDevEigen;     // maximum absolute deviation of diagonal elements from 1
  int64 tix;
  SPCCalculateCallbackData spccb; // passed to caller
  double* CovMatrix;          // mutual covariation           N*N | M*M
  double* EigenVec;           // eigenvectors of matrix       N*N | M*M
  double* EigenVal;           // eigen values of matrix       N   | M
  pfnCovariation pfnCov;
  void* extbuftmp,*bup = NULL;  // temporary buffer
  int nTMP;                     // size of temporary buffer

    // check arguments
    if (pSPCB==NULL)
      return ERR_GEN_NULLPOINTER;
    for (i=0;i<pSPCB->nSampleCount;i++)
      if (!pSPCB->Samples_Present[i])
        return MIA_PCA_NOSAMPLES;
    // make fake if necessary
    if (pCB==NULL)
      pCB = ownBPL_fnFakePCCalculateCallback;
    // pass callback token to callback structure
    spccb.pvToken = pvToken;
    // choosing a per-sample or per-element covariation
    if (pSPCB->isPerElem)
    { // per-element covariation
      eigendim = pSPCB->nSampleDimension;
      pfnCov = &ownSPCB_CovariatePerElements;
      strcpy(spccb.repstr,"PCCalc: per-element covariation mode");
    }
    else
    { // per-sample covariation
      eigendim = pSPCB->nSampleCount;
      pfnCov = &ownSPCB_CovariatePerSamples;
      strcpy(spccb.repstr,"PCCalc: per-sample covariation mode");
    }
    (*pCB)(&spccb);
    // allocate memory
    if ((bup = malloc(
           sizeof(CovMatrix[0])*eigendim*eigendim+
           sizeof(EigenVec[0] )*eigendim*eigendim+
           sizeof(EigenVal[0] )*eigendim+
           (nTMP =
           sizeof(int)*pSPCB->nSampleDimension+
           sizeof(double)*pSPCB->nSampleCount)))==NULL)
      return ERR_GEN_NOMEMORY;
    CovMatrix = (double*)bup;
    EigenVec = CovMatrix+eigendim*eigendim;
    EigenVal = EigenVec+eigendim*eigendim;
    extbuftmp = EigenVal+eigendim;
    // calculate mean vector and normalise image set to a zero mean
    strcpy(spccb.repstr,"PCCalc: running NormaliseOnMean");
    (*pCB)(&spccb);
    tix = clock();
    if ((res = ownSPCB_NormaliseOnMean(
      pSPCB->MeanSample,
      pSPCB->Samples,
      pSPCB->nSampleDimension,
      pSPCB->nSampleCount,
      extbuftmp,
      nTMP))!=ERR_OK)
      goto endpcbcalc;
    tix = clock()-tix;
    sprintf(spccb.repstr,"  Done in %d clox",(int)tix);
    (*pCB)(&spccb);
    // calculate per-element or per-sample covariation
    strcpy(spccb.repstr,"PCCalc: Calculating covariation matrix");
    (*pCB)(&spccb);
    tix = clock();
    (*pfnCov)(
      CovMatrix,
      pSPCB->Samples,
      pSPCB->nSampleDimension,
      pSPCB->nSampleCount);
    tix = clock()-tix;
    memcpy(EigenVec,CovMatrix,eigendim*eigendim*sizeof(EigenVec[0]));
    sprintf(spccb.repstr,"  Done in %d clox",(int)tix);
    (*pCB)(&spccb);
    // eigen vector calculation
    strcpy(spccb.repstr,"PCCalc: running QL eigen");
    (*pCB)(&spccb);
    tix = clock();
    if ((res = BPL_Matrix_Eigen_Symm(
      EigenVec, // IN:  source matrix / OUT: eigenvector matrix
      EigenVal, // OUT: eigen values
      eigendim,  // IN:  size of matrix
      extbuftmp,  // IN:  temporary buffer
      nTMP))!=ERR_OK)
      goto endpcbcalc;
    tix = clock()-tix;
    sprintf(spccb.repstr,"  Done in %d clox",(int)tix);
    (*pCB)(&spccb);
    // vectors are in eigen-increase order. reverse the order
    for (i=0;i<eigendim/2;i++)
    {
      // exchange vectors
      memcpy(extbuftmp,
             EigenVec+eigendim*i,
             eigendim*sizeof(EigenVec[0]));
      memcpy(EigenVec+eigendim*i,
             EigenVec+eigendim*(eigendim-i-1),
             eigendim*sizeof(EigenVec[0]));
      memcpy(EigenVec+eigendim*(eigendim-i-1),
             extbuftmp,
             eigendim*sizeof(EigenVec[0]));
      // exchange values
      ErrMaxAbsDevDiag = EigenVal[i];
      EigenVal[i] = EigenVal[eigendim-i-1];
      EigenVal[eigendim-i-1] = ErrMaxAbsDevDiag;
    }
    // calc eigen errors
    strcpy(spccb.repstr,"PCCalc: checking eigen errors");
    (*pCB)(&spccb);
    BPL_Matrix_EstimateErrorNonUnitar(EigenVec,eigendim,
      &ErrMaxAbsDevDiag,&ErrMaxAbsDevNdiag,
      &ErrAvgAbsDevDiag,&ErrAvgAbsDevNdiag);
    sprintf(spccb.repstr,"  Eigen matrix non-unitarity:\n\tMaxDiag = %le,\n\tMaxNdiag = %le,\n\tAvgDiag = %le,\n\tAvgNdiag = %le",
      ErrMaxAbsDevDiag,ErrMaxAbsDevNdiag,ErrAvgAbsDevDiag,ErrAvgAbsDevNdiag);
    (*pCB)(&spccb);
    BPL_Matrix_EstimateErrorNonEigen(
      CovMatrix,    // IN:  source matrix
      EigenVec,  // IN:  eigen vectors
      EigenVal,  // IN:  eigen values
      eigendim,  // IN:  size of source matrix
      pSPCB->nNecessaryVectors,  // IN:  number of eigen vectors to be tested
      &ErrMaxDevEigen);
    sprintf(spccb.repstr,"  Maximum of ||Ax-Lx||/L = %le",ErrMaxDevEigen);
    (*pCB)(&spccb);
    // quit if big error
    if ((ErrMaxAbsDevDiag*pSPCB->nSampleDimension*pSPCB->nSampleCount>1.)||
        (ErrMaxAbsDevNdiag*pSPCB->nSampleDimension*pSPCB->nSampleCount>1.)||
        (ErrMaxDevEigen*pSPCB->nSampleDimension*pSPCB->nSampleCount>1.))
    {
      strcpy(spccb.repstr,"WARNING! Eigen discrepancy too big.");
      (*pCB)(&spccb);
//      return MIA_PCA_TOO_BIG_DISCREPANCY;
    }
    // find the order of first eigen value that is less than M/sqrt(12)
    ErrMaxAbsDevDiag = eigendim/sqrt(12.);
    for (i=0;i<eigendim;i++)
      if (EigenVal[i]<=ErrMaxAbsDevDiag)
        break;
    // report on eigen vector border
    if (i==eigendim)
      // this is very unlikely to happen
      strcpy(spccb.repstr,"All eigen values exceed M/sqrt(12)");
    else
      sprintf(spccb.repstr,"L[%d] = %f < %f = M/sqrt(12)",
        i,(float)(EigenVal[i]),(float)ErrMaxAbsDevDiag);
    (*pCB)(&spccb);
    // quit if big error
    if (i<=pSPCB->nNecessaryVectors)
    {
      strcpy(spccb.repstr,"WARNING! PC space vectors may be strongly influenced by rounding error.");
      (*pCB)(&spccb);
//      return MIA_PCA_INSUFFICIENT_RELIABLE_EIGENS_ROUNDING;
    }
    // small and negative eigen values filtration
    startclear = eigendim;
    for (i=1;i<eigendim;i++)
      if (EigenVal[0]/EigenVal[i]>eigendim*eigendim)
      {
        startclear = i;
        break;
      }
    // maybe negative or zero can cause error. make positive
    for (i=startclear;i<eigendim;i++)
    {
      memset(EigenVec+i*eigendim,0,eigendim*sizeof(double));
      EigenVal[i] = EigenVal[startclear-1];
    }
    sprintf(spccb.repstr,"PCCalc: last %d of %d eigens and vectors are zeroed",eigendim-startclear,eigendim);
    (*pCB)(&spccb);
    // quit if big error
    if (startclear<=pSPCB->nNecessaryVectors)
    {
      strcpy(spccb.repstr,"WARNING! Eigen vectors are poorly conditioned.");
      (*pCB)(&spccb);
//      return MIA_PCA_INSUFFICIENT_RELIABLE_EIGENS_CONDITIONALITY;
    }
    // obtain the space
    if (pSPCB->isPerElem)
    {
      memcpy(pSPCB->EigenSpace,EigenVec,
        pSPCB->nNecessaryVectors*pSPCB->nSampleDimension*sizeof(pSPCB->EigenSpace[0]));
    }
    else
    {
      strcpy(spccb.repstr,"PCCalc: running MatrixSpecProductABL");
      (*pCB)(&spccb);
      tix = clock();
      ownSPCB_MatrixSpecProductABL(
        pSPCB->EigenSpace,         // destination: double matrix IxK
        EigenVec,   // source1: double matrix     IxJ
        pSPCB->Samples,     // source2: char matrix       JxK
        EigenVal,    // source3: values            I
        pSPCB->nNecessaryVectors, // size I
        pSPCB->nSampleCount,      // size J
        pSPCB->nSampleDimension); // size K
      tix = clock()-tix;
      sprintf(spccb.repstr,"  Done in %d clox",(int)tix);
      (*pCB)(&spccb);
    }
    // ok
    strcpy(spccb.repstr,"PCCalc: finished");
    (*pCB)(&spccb);
endpcbcalc:
    if (bup)
      free(bup);
    return res;
}

// size of buffer is: PDim
MIA_RESULT_CODE BPL_PCB_SavePC(
  const SPCBuilder* pSPCB,      // IN: PC space to be saved
  unsigned char* dst,           // IN: buffer accepting space
  int* dstsize,                 // IN: size of buffer, OUT: bytes actually placed in buffer
  int PDim)                     // IN: space dimension
{
  int IDim,i,j;
  double maxmod;
  unsigned char *shs = NULL,*orig_dst = dst;
  MIA_RESULT_CODE res = ERR_OK;

    // check arguments
    if ((pSPCB==NULL)||(dst==NULL)||(dstsize==NULL))
      return ERR_GEN_NULLPOINTER;
    if (PDim%PC_NVEC_ALIGN)
      return ERR_GEN_BAD_ALIGNMENT;
    if (PDim>pSPCB->nNecessaryVectors)
      return ERR_GEN_INVALID_PARAMETER;
    // calc IDim to be a multipe of SUMCHUNK
    IDim = ((pSPCB->nSampleDimension+SUMCHUNK-1)/SUMCHUNK)*SUMCHUNK;
    // ID
    *((unsigned int*)dst) = MEMID_SPCANALYSER;
    dst += 4;
    // dimensions
    *((short*)dst) = (short)IDim;
    dst += 2;
    *((short*)dst) = (short)PDim;
    dst += 2;
    // mean sample
    memcpy(dst,pSPCB->MeanSample,pSPCB->nSampleDimension);
    dst += pSPCB->nSampleDimension;
    for (i=IDim-pSPCB->nSampleDimension;i;i--)
      *dst++ = 0;
    // eigen space
    if ((shs = (unsigned char*)malloc(PDim))==NULL)
      goto endpcbsave;
    // find shifts
    for (i=0;i<PDim;i++)
    {
      maxmod = 0.;
      for (j=0;j<IDim;j++)
      {
        if (pSPCB->EigenSpace[i*IDim+j]>maxmod)
          maxmod = pSPCB->EigenSpace[i*IDim+j];
        else
          if (-pSPCB->EigenSpace[i*IDim+j]>maxmod)
            maxmod = -pSPCB->EigenSpace[i*IDim+j];
      }
      if (maxmod==0.)
      {
        res = ERR_GEN_INTERNAL;  // eigenspace vector is wholy zero
        goto endpcbsave;
      }
      shs[i] = (unsigned char)MIA_IntLog2Hi((unsigned int)((1<<(SHORTSIZE-1))/maxmod));
    }
    // write vectors
    for (i=0;i<pSPCB->nSampleDimension;i++)
      for (j=0;j<PDim;j++)
      {
        *((short*)dst) = (short)(pSPCB->EigenSpace[j*IDim+i]*(1<<shs[j]));
        dst += 2;
      }
    // zero space up to SUMCHUNK*PDim alignment
    for (i=(IDim-pSPCB->nSampleDimension)*PDim;i;i--)
    {
      *((short*)dst) = 0;
      dst += 2;
    }
    // write shifts
    for (j=0;j<PDim;j++)
      *dst++ = shs[j];
    // output number of bytes saved
    *dstsize = (int)(dst-orig_dst);
endpcbsave:
    if (shs)
      free(shs);
    return res;
}

MIA_RESULT_CODE BPL_PCB_Expand(
  const SPCBuilder* pSPCB,  // IN:  PCA structure
  double* coefs,            // OUT: expansion coefficients
  const unsigned char* im)  // IN:  source image
{
  int i,j;
  double sum;

    for (i=0;i<pSPCB->nNecessaryVectors;i++)
    {
      sum = 0;
      for (j=0;j<pSPCB->nSampleDimension;j++)
        sum += pSPCB->EigenSpace[pSPCB->nSampleDimension*i+j]*(((int)((unsigned short)(im[j])))-((int)((unsigned short)(pSPCB->MeanSample[j]))));
      coefs[i] = sum;
    }
    return ERR_OK;
}
