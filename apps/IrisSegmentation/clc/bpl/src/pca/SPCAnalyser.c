/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:     21 March 2001                                       */
/*      Modified:    5 November 2001 - small (and negative) eigen values */
/*                                     filtration added                  */
/*      Modified:    6 November 2001 - added cutted save/load (v2)       */
/*      Modified:    11 May 2004 - ported to C, double->int (v3)         */
/*      Revision:    3.0.00                                              */
/*      Purpose:     principal components analyser                       */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/
#define __FILENUM__ 98 // __FILENUM__TAG98

#include "bpl.h"

// set structure according to block
MIA_RESULT_CODE BPL_PCA_ArrangePCStruct(
  SPCAnalyser* pSPCA,       // IN:  structure to be set
  const void* inbuf,        // IN:  buffer storing the structure
  int* nBytes)              // OUT: number of bytes used in buffer
{
  unsigned char *ib;

    // check arguments
    if ((pSPCA==NULL)||(inbuf==NULL))
      return ERR_GEN_NULLPOINTER;
    // check ID
    ib = (unsigned char*)inbuf;
    if (*((unsigned int*)ib)!=MEMID_SPCANALYSER)
      return ERR_GEN_INVALID_OBJECT_ID;
    ib += 4;    // p+4
    // set dimensions
    pSPCA->IDim = *((short*)ib);
    ib += 2;    // p+6
    pSPCA->PDim = *((short*)ib);
    ib += 2;    // p+8
    // set data
    pSPCA->Mean = ib;
    ib += pSPCA->IDim;
    pSPCA->Vecs = (short*)ib;
    ib += sizeof(short)*pSPCA->IDim*pSPCA->PDim;
    pSPCA->Shifts = ib;
    *nBytes = 4+2+2+                                  // header
              pSPCA->IDim+                            // mean
              sizeof(short)*pSPCA->IDim*pSPCA->PDim+  // vecs
              pSPCA->PDim;                            // shifts
    return ERR_OK;
}

// expand the incoming data
MIA_RESULT_CODE BPL_PCA_Expand(
  const SPCAnalyser* pSPCA, // IN:  PCA structure
  short* coefs,             // OUT: expansion coefficients
  const unsigned char* im,  // IN:  source image
  void* extbuftmp)
{
/*  simple analog
  short v;
  const short *vex;
  const unsigned char* mean;
  int i,k,pd;
  int64* Zums;

    // init
    mean = pSPCA->Mean;
    pd = pSPCA->PDim;
    Zums = (int64*)extbuftmp;
    // second part of sums is for totals
    for (k=0;k<pd;k++)
      Zums[k] = 0;
    // loop
    for (i=0;i<pSPCA->IDim;i++)
    {
      v = (short)(((short)((unsigned short)(im[i])))-((short)((unsigned short)(mean[i]))));
      vex = pSPCA->Vecs+i*pd;
      // loop by vector components
      for (k=0;k<pd;k++)
        Zums[k] += v*vex[k];
    }
    // draw coefs from total sums
    for (k=0;k<pd;k++)
      coefs[k] = (short)((Zums[k])>>(pSPCA->Shifts[k]));
    return ERR_OK;
*/
  short v;
  const short *vex;
  const unsigned char* mean;
  int i,j,k,pd;
  int *sums;

    // check arguments
    if ((pSPCA==NULL)||(coefs==NULL)||(im==NULL)||(extbuftmp==NULL))
      return ERR_GEN_NULLPOINTER;
    // init
    mean = pSPCA->Mean;
    pd = pSPCA->PDim;
    sums = (int*)extbuftmp;
    // second part of sums is for totals
    for (k=pd-1;k>=0;k--)
      sums[pd+k] = 0;
    // loop by chunks
    for (i=pSPCA->IDim-SUMCHUNK;i>=0;i-=SUMCHUNK)
    {
      // clear chunk sums
      for (k=pd-1;k>=0;k--)
        sums[k] = 0;
      // loop inside chunk
      for (j=i+SUMCHUNK-1;j>=i;j--)
      {
        v = (short)(((short)((unsigned short)(im[j])))-((short)((unsigned short)(mean[j]))));
        vex = pSPCA->Vecs+j*pd;
        // loop by vector components
        for (k=pd-1;k>=0;k--)
          sums[k] += v*vex[k];
      }
      // add chunk sums to total sums
      for (k=pd-1;k>=0;k--)
        if (sums[k]>0)
          sums[pd+k] += sums[k]>>CHARSIZE;
        else
          sums[pd+k] += (int)((((unsigned int)(sums[k]))>>CHARSIZE)+(0u-(1<<(INTSIZE-CHARSIZE))));
    }
    // draw coefs from total sums
    for (k=pd-1;k>=0;k--)
      coefs[k] = (short)((sums[pd+k])>>(pSPCA->Shifts[k]-CHARSIZE));
    return ERR_OK;
}
