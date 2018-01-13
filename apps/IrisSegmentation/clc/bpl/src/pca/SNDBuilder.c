/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:     5 December 2007                                     */
/*      Revision:    1.0.00                                              */
/*      Purpose:     Normal distribution classifier builder              */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#define __FILENUM__ 44 // __FILENUM__TAG44

#include <math.h>
#include <malloc.h>
#include <memory.h>
#include "SNDClassifier.h"
#include "bpl.h"

typedef struct
{
  int nClasses;
  int nDimensions;
  int* pnSamplesInClasses;
  int* pnSamplesAllocated;
  double** pdSamples;
} SNDClassBuilder;

// initialise classifier builder
MIA_RESULT_CODE BPL_NDC_Create(
  HNDClassBuilder* phNDB, // OUT: handle to classifier builder
  int nClasses,           // IN:  number of classes
  int nDimensions)        // IN:  number of dimensions
{
  int sz;
  SNDClassBuilder* psNDB;

    // check arguments
    if (phNDB==NULL)
      return ERR_GEN_NULLPOINTER;
    if ((nClasses<=1)||(nDimensions<=1)||(nClasses>0xffff)||(nDimensions>0xffff))
      return ERR_GEN_INVALID_PARAMETER;
    // allocate
    if ((psNDB = (SNDClassBuilder*)malloc(sz = 
                    sizeof(psNDB[0])+
                    nClasses*sizeof(psNDB->pnSamplesInClasses[0])+
                    nClasses*sizeof(psNDB->pnSamplesAllocated[0])+
                    nClasses*sizeof(psNDB->pdSamples[0])))==NULL)
      return ERR_GEN_NOMEMORY;
    memset(psNDB,0,sz);
    psNDB->pnSamplesInClasses = (int*)(psNDB+1);
    psNDB->pnSamplesAllocated = psNDB->pnSamplesInClasses+nClasses;
    psNDB->pdSamples = (double**)(psNDB->pnSamplesAllocated+nClasses);
    // set vars
    psNDB->nClasses = nClasses;
    psNDB->nDimensions = nDimensions;
    // output
    *phNDB = (HNDClassBuilder)psNDB;
    return ERR_OK;
}

// free builder instance
MIA_RESULT_CODE BPL_NDC_FreeBuilder(
  HNDClassBuilder* phNDB) // IN:  handle to classifier
{
  int i;
  SNDClassBuilder* psNDB;

    if (phNDB==NULL)
      return ERR_GEN_NULLPOINTER;
    if (*phNDB==NULL)
      return ERR_GEN_NULLPOINTER;
    psNDB = (SNDClassBuilder*)(*phNDB);
    for (i=0;i<psNDB->nClasses;i++)
      if (psNDB->pdSamples[i])
        free(psNDB->pdSamples[i]);
    free(*phNDB);
    *phNDB = NULL;
    return ERR_OK;
}

// add sample to teaching set
MIA_RESULT_CODE BPL_NDC_AddSample(
  HNDClassBuilder hNDB,   // IN:  classifier builder handle
  const double* pdVec,    // IN:  sample
  int nClass)             // IN:  class number
{
  SNDClassBuilder* psNDB = (SNDClassBuilder*)hNDB;

    // check arguments
    if ((hNDB==NULL)||(pdVec==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((nClass<0)||(nClass>=psNDB->nClasses))
      return ERR_GEN_INVALID_PARAMETER;
    // do we need more memory?
    if (psNDB->pnSamplesInClasses[nClass]==psNDB->pnSamplesAllocated[nClass])
      // realloc memory
      if ((psNDB->pdSamples[nClass] = (double*)realloc(
              psNDB->pdSamples[nClass],
              (psNDB->pnSamplesAllocated[nClass] += 1024)*psNDB->nDimensions*
                sizeof(psNDB->pdSamples[nClass][0])))==NULL)
        return ERR_GEN_NOMEMORY;
    // copy data
    memcpy(psNDB->pdSamples[nClass]+
             psNDB->pnSamplesInClasses[nClass]*psNDB->nDimensions,
           pdVec,psNDB->nDimensions*sizeof(pdVec[0]));
    psNDB->pnSamplesInClasses[nClass]++;
    return ERR_OK;
}

// build the classifier
MIA_RESULT_CODE BPL_NDC_Build(
  HNDClassBuilder hNDB,   // IN:  classifier builder handle
  void* pvBuf,            // OUT: buffer for classifier packed data
  int* pnLen,             // IN/OUT: bytes allocated/used for storing
  const double* pdPens)   // IN:  penalty matrix
{
  int i,sz,j,k,l,dim,cla;
  SNDClassBuilder* psNDB = (SNDClassBuilder*)hNDB;
  double *means,*matrices,*invmatrices,*packmatrices,*S,*m;
  unsigned char* pucbuf;

    // check arguments
    if ((hNDB==NULL)||(pnLen==NULL))
      return ERR_GEN_NULLPOINTER;
    for (i=0;i<psNDB->nClasses;i++)
      if (psNDB->pnSamplesInClasses[i]<=2)
        return ERR_GEN_NOT_INITIALISED;
    if (((ptr_t)pvBuf)&3)
      return ERR_GEN_BAD_ALIGNMENT;
    // handy variables
    dim = psNDB->nDimensions;
    cla = psNDB->nClasses;
    // calc size
    sz =  sizeof(int)+                      // ID
          2*sizeof(short)+                  // sizes
          sizeof(float)*cla*dim+            // means
          sizeof(float)*cla*dim*(dim+1)/2+  // packed matrices
          sizeof(float)*cla*(cla-1);        // penalties
    if (pvBuf==NULL)
    {
      *pnLen = sz;
      return ERR_OK;
    }
    if (*pnLen<sz)
    {
      *pnLen = sz;
      return ERR_GEN_INSUFFICIENT_BUFFER;
    }
    // allocate means and matrices
    means = (double*)malloc(sizeof(means[0])*cla*
              (dim+3*dim*dim));
    if (means==NULL)
      return ERR_GEN_NOMEMORY;
    memset(means,0,sizeof(means[0])*cla*
              (dim+3*dim*dim));
    matrices = means+cla*dim;
    packmatrices = matrices+cla*dim*dim;
    invmatrices = packmatrices+cla*dim*dim;
    // calculate means and matrices
    for (i=0;i<cla;i++)                 // for each class
    {
      m = means+i*dim;
      S = matrices+i*dim*dim;
      // sums
      for (j=0;j<psNDB->pnSamplesInClasses[i];j++)  // for each sample
        for (k=0;k<dim;k++)          // for each dimension
          m[k] += psNDB->pdSamples[i][j*dim+k];
      // means
      for (j=0;j<dim;j++)          // for each dimension
        m[j] /= psNDB->pnSamplesInClasses[i];
      // matrices. some excess calculations here. can omit almost half
      // of calcs since matrix is symmetric
      for (j=0;j<dim;j++)            // for each dimension
        for (k=0;k<dim;k++)          // for each dimension
          for (l=0;l<psNDB->pnSamplesInClasses[i];l++)
            S[j*dim+k] += 
              psNDB->pdSamples[i][l*dim+j]*
              psNDB->pdSamples[i][l*dim+k];
      // normalize
      for (j=0;j<dim;j++)            // for each dimension
        for (k=0;k<dim;k++)          // for each dimension
          S[j*dim+k] = S[j*dim+k]/psNDB->pnSamplesInClasses[i]-m[j]*m[k];
    }
    // pack matrices
    for (i=0;i<cla;i++)                 // for each class
      BPL_Matrix_Pack(
        packmatrices+i*dim*dim,
        matrices+i*dim*dim,
        dim);
    // invert
    for (i=0;i<cla;i++)                 // for each class
      if (BPL_Matrix_InvertPacked(
            invmatrices+i*dim*dim,
            packmatrices+i*dim*dim,
            dim))
      {
        free(means);
        return ERR_GEN_SIZE_NOT_MATCH;
      }
/*    // unpack
    for (i=0;i<cla;i++)                 // for each class
      BPL_Matrix_Pack(
        matrices+i*dim*dim,
        invmatrices+i*dim*dim,
        dim);*/
    // pack now!
    pucbuf = (unsigned char*)pvBuf;
    // ID
    *((int*)pucbuf) = MEMID_NORMDISTRCLASS;
    pucbuf += sizeof(int);
    // class number
    *((short*)pucbuf) = (short)(cla);
    pucbuf += sizeof(short);
    // dimension number
    *((short*)pucbuf) = (short)(dim);
    pucbuf += sizeof(short);
    // means
    for (i=0;i<cla*dim;i++)
      ((float*)pucbuf)[i] = (float)(means[i]);
    pucbuf += sizeof(float)*cla*dim;
    // packed matrices
    for (i=0;i<cla;i++)
    {
      for (j=0;j<dim*(dim+1)/2;j++)
        ((float*)pucbuf)[j] = (float)(invmatrices[i*dim*dim+j]);
      pucbuf += sizeof(float)*dim*(dim+1)/2;
    }
    // penalties
    if (pdPens)
    {
      for (i=0;i<cla;i++)
        for (j=0;j<cla;j++)
          if (i!=j)
          {
            *((float*)pucbuf) = (float)(pdPens[i*cla+j]);
            pucbuf += sizeof(float);
          }
    }
    else
    {
      for (i=0;i<cla*(cla-1);i++)
        ((float*)pucbuf)[i] = 1.f;
      pucbuf += sizeof(float)*cla*(cla-1);
    }
    // clear memstock
    free(means);
    if (sz!=(int)((ptr_t)pucbuf-(ptr_t)pvBuf))
      return ERR_GEN_INTERNAL;
    *pnLen = sz;
    return ERR_OK;
}
