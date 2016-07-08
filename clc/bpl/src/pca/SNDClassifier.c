/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:     5 December 2007                                     */
/*      Revision:    1.0.00                                              */
/*      Purpose:     Normal distribution classifier builder              */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#include <math.h>
#include "SNDClassifier.h"
#include "bpl.h"

typedef struct
{
  void* pvMembuf;
  int nClasses;
  int nDimensions;
  float* pdMeans;   // means of classes
  float* pdCorrs;   // inverted covariations
  float* pdPenls;   // penalties
  float* Fj;        // temporary array of nclasses elements
} SNDClassifier;

// initialise classifier structure from data set by ndc_build function
MIA_RESULT_CODE BPL_NDC_Initialise(
  HNDClassifier* phNDC, // OUT: handle to classifier
  int* pnValNumber,     // OUT: number of variables
  const void* pvPacked, // IN:  packed classifier data (can be cleaned after this function exits)
  void* pvMembuf,       // IN:  memory buffer for unpacking classifier
  int* pnBuflen)        // IN/OUT: number of bytes allocated/used in buffer
{
  int sz,dim,cla,i,j;
  unsigned char* pucbuf;
  SNDClassifier* pSC;

    // check arguments
    if ((pvPacked==NULL)||(pnBuflen==NULL))
      return ERR_GEN_NULLPOINTER;
    if (((ptr_t)pvMembuf)&3)
      return ERR_GEN_BAD_ALIGNMENT;
    // unpack dimensions
    pucbuf = (unsigned char*)pvPacked;
    // ID
    if (*((int*)pucbuf)!=MEMID_NORMDISTRCLASS)
      return ERR_GEN_INVALID_OBJECT_ID;
    pucbuf += sizeof(int);
    // class number
    cla = *((short*)pucbuf);
    pucbuf += sizeof(short);
    // dimension number
    dim = *((short*)pucbuf);
    pucbuf += sizeof(short);
    // check validity
    if ((cla<=1)||(dim<=1))
      return ERR_GEN_INVALID_PARAMETER;
    // calc size
    sz = sizeof(SNDClassifier)+
         sizeof(pSC->pdMeans[0])*cla*dim+
         sizeof(pSC->pdCorrs[0])*cla*dim*dim+
         sizeof(pSC->pdPenls[0])*cla*cla+
         sizeof(pSC->Fj[0])     *cla;
    if (pvMembuf==NULL)
    {
      *pnBuflen = sz;
      return ERR_OK;
    }
    if (*pnBuflen<sz)
    {
      *pnBuflen = sz;
      return ERR_GEN_INSUFFICIENT_BUFFER;
    }
    *pnBuflen = sz;
    // additional check
    if (phNDC==NULL)
      return ERR_GEN_NULLPOINTER;
    // allocate
    pSC = (SNDClassifier*)pvMembuf;
    pSC->pdMeans = (float*)(pSC+1);
    pSC->pdCorrs = pSC->pdMeans+cla*dim;
    pSC->pdPenls = pSC->pdCorrs+cla*dim*dim;
    pSC->Fj      = pSC->pdPenls+cla*cla;
    // means
    for (i=0;i<cla*dim;i++)
      pSC->pdMeans[i] = ((float*)pucbuf)[i];
    pucbuf += sizeof(float)*cla*dim;
    // packed matrices
    for (i=0;i<cla;i++)
    {
      for (j=0;j<dim*(dim+1)/2;j++)
        pSC->pdCorrs[i*dim*dim+j] = ((float*)pucbuf)[j];
      pucbuf += sizeof(float)*dim*(dim+1)/2;
      // unpack
      BPL_Matrix_UnpackFloat(pSC->pdCorrs+i*dim*dim,pSC->pdCorrs+i*dim*dim,dim);
    }
    // penalties
    for (i=0;i<cla;i++)
      for (j=0;j<cla;j++)
        if (i!=j)
        {
          pSC->pdPenls[i*cla+j] = *((float*)pucbuf);
          pucbuf += sizeof(float);
        }
        else
          pSC->pdPenls[i*cla+j] = 0.f;
    // set vars
    pSC->nClasses = cla;
    pSC->nDimensions = dim;
    pSC->pvMembuf = pvMembuf;
    // output
    *phNDC = (HNDClassifier)pSC;
    if (pnValNumber)
     *pnValNumber = dim;
    return ERR_OK;
}

// free classifier instance
MIA_RESULT_CODE BPL_NDC_Free(
  void** ppvMembuf,     // OUT: membuf used for classifier
  HNDClassifier* phNDC) // IN:  handle to classifier
{
  SNDClassifier* pSC;

    if (phNDC==NULL)
      return ERR_GEN_NULLPOINTER;
    pSC = (SNDClassifier*)(*phNDC);
    if (pSC==NULL)
      return ERR_GEN_NOT_INITIALISED;
    *phNDC = NULL;
    if (ppvMembuf)
      *ppvMembuf = pSC->pvMembuf;
    return ERR_OK;
}

#define NUMINOZO

// classify the sample
MIA_RESULT_CODE BPL_NDC_Classify(
  int* pnClass,         // OUT: class number
  float* pdRisks,       // OUT: risk value for each class
  HNDClassifier hNDC,   // IN:  classifier handle
  const float* pdVec_)  // OUT: vector to be classified
{
  int i,j,k,cla,dim,minidx;
  SNDClassifier* pSC;
  float F,S,*C,*m,minrisk;

float pdVec[32];

    // check arguments
    if ((hNDC==NULL)||(pdVec_==NULL)||((pnClass==NULL)&&(pdRisks==NULL)))
      return ERR_GEN_NULLPOINTER;
    pSC = (SNDClassifier*)hNDC;
    // handy vars
    cla = pSC->nClasses;
    dim = pSC->nDimensions;
    MIA_memcpy(pdVec,pdVec_,sizeof(pdVec[0])*dim);
#ifdef NUMINOZO
    if (cla==2)
    {
      j = k = 0;
      for (i=0;i<dim;i++)
      {
        if (pSC->pdMeans[i]<pSC->pdMeans[i+dim])
        {
          if (pdVec[i]<pSC->pdMeans[i])
            j++;
          else
            if (pSC->pdMeans[i+dim]<pdVec[i])
              k++;
        }
        else
          if (pdVec[i]>pSC->pdMeans[i])
            k++;
          else
            if (pSC->pdMeans[i+dim]>pdVec[i])
              j++;
      }
      if ((j+k>3)&&(!(j*k)))
      {
        for (i=0;i<dim;i++)
        {
          if (pSC->pdMeans[i]<pSC->pdMeans[i+dim])
          {
            if (pdVec[i]<pSC->pdMeans[i])
              pdVec[i] = pSC->pdMeans[i];
            else
              if (pSC->pdMeans[i+dim]<pdVec[i])
                pdVec[i] = pSC->pdMeans[i+dim];
          }
          else
            if (pdVec[i]>pSC->pdMeans[i])
              pdVec[i] = pSC->pdMeans[i];
            else
              if (pSC->pdMeans[i+dim]>pdVec[i])
                pdVec[i] = pSC->pdMeans[i+dim];
        }
      }
    }
#endif //NUMINOZO
    // calculate class as: K = arg(min(Ri))
    // Ri = Sum(Fj*rij)   / by all j!=i
    // Fj = (x-mj)' * Cj * (x-mj) / quadratic form, x,mj - vectors, Cj - matrix
    // for each class
    for (i=0;i<cla;i++)
    {
      m = pSC->pdMeans+i*dim;
      C = pSC->pdCorrs+i*dim*dim;
      F = 0;
      for (j=0;j<dim;j++)
      {
        S = 0;
        for (k=0;k<dim;k++)
          S += C[j*dim+k]*(pdVec[k]-m[k]);
        F += S*(pdVec[j]-m[j]);
      }
      pSC->Fj[i] = F;
    }
    // look for minimal risk
    minidx = -1;
    minrisk = 1e10;
    for (i=0;i<cla;i++)
    {
      F = 0;
      for (j=0;j<cla;j++)
//        F -= pSC->Fj[j]*pSC->pdPenls[i*cla+j];
        F += (float)(exp(-pSC->Fj[j])*pSC->pdPenls[i*cla+j]);
//        F += (float)(pSC->pdPenls[i*cla+j]/(1.+pSC->Fj[j]*pSC->Fj[j]));
      if (pdRisks)
        pdRisks[i] = F;
      if (minidx<0)
      {
        minrisk = F;
        minidx = 0;
      }
      else
        if (minrisk>F)
        {
          minrisk = F;
          minidx = i;
        }
    }
    if (pnClass)
      *pnClass = minidx;
    return ERR_OK;
}
