/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  19 March 2002                                          */
/*      Revision: 1.0.00                                                 */
/*      Purpose:  Fast Walsh and Haar transforms                         */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#ifdef DATATYPE
#include <memory.h>
#include <malloc.h>
#include "ipl.h"

// perform one butterfly level of Haar transform
staticFunc void ownIPL_FHT_Stage(
  DATATYPE* dst,
  const DATATYPE* src,
  unsigned int halfsize)
{
  unsigned int i = halfsize;

    while (i--)
    {
      dst[i         ] = src[2*i]+src[2*i+1];
      dst[i+halfsize] = src[2*i]-src[2*i+1];
    }
}

// transforming 1-dimensional data with power of 2 length
staticFunc void ownIPL_FastHaarTransform1D(
  DATATYPE* dst,        // destination
  DATATYPE* tmp,        // temporary array
  unsigned int lenpow) // log2(length)
{
  int i;
  DATATYPE *t;

    // butterflies
    for (i=1<<(lenpow-1);i;i>>=1)
    // this cycle runs 'lenpow' times
    {
      ownIPL_FHT_Stage(tmp,dst,i);
      t   = dst;
      dst = tmp;
      tmp = t;
    }
    // collect data to one place
    if (lenpow&1)
    { // tmp and dst were changed
      for (i=lenpow-1;i>=0;i-=2)
        memcpy(tmp+(((size_t)1)<<i),dst+(((size_t)1)<<i),(((size_t)1)<<i)*sizeof(DATATYPE));
      *tmp = *dst;
    }
    else
      // tmp and dst are original
      for (i=lenpow-1;i>0;i-=2)
        memcpy(dst+(((size_t)1)<<i),tmp+(((size_t)1)<<i),(((size_t)1)<<i)*sizeof(DATATYPE));
}

// transforming 2-dimensional data with power of 2 sizes
staticFunc void ownIPL_FastHaarTransform2D(
  DATATYPE* dst,          // destination
  DATATYPE* tmp,          // temporary array
  unsigned int lenpow_x, // length along x
  unsigned int lenpow_y) // length along y
{
  unsigned int i,j,size_x,size_y;

    size_x = 1<<lenpow_x;
    size_y = 1<<lenpow_y;
    // x transform
    for (i=0;i<size_y;i++)
      ownIPL_FastHaarTransform1D(dst+i*size_x,tmp,lenpow_x);
    // flip image
    for (i=0;i<size_y;i++)
      for (j=0;j<size_x;j++)
        tmp[i+j*size_y] = dst[i*size_x+j];
    // y transform
    for (j=0;j<size_x;j++)
      ownIPL_FastHaarTransform1D(tmp+j*size_y,dst,lenpow_y);
    // flip image back
    for (i=0;i<size_y;i++)
      for (j=0;j<size_x;j++)
        dst[i*size_x+j] = tmp[i+j*size_y];
}

// 1D Haar transform of doubles
MIA_RESULT_CODE IPL_FastHaarTransform1D(
  DATATYPE* dst,          // destination
  DATATYPE* tmp,          // temporary. may be NULL
  unsigned int lenpow)   // log2(length)
{
  double *t;

    if (dst==NULL)
      return ERR_GEN_NULLPOINTER;
    if ((lenpow<2)||(lenpow>16))
      return ERR_GEN_INVALID_PARAMETER;
    if (tmp!=NULL)
      t = tmp;
    else
      if ((t = (double*)malloc((((size_t)1)<<lenpow)*sizeof(double)))==NULL)
        return ERR_GEN_NOMEMORY;
    ownIPL_FastHaarTransform1D(dst,t,lenpow);
    if (tmp==NULL)
      free(t);
    return ERR_OK;
}

// 2D Haar transform of doubles
MIA_RESULT_CODE IPL_FastHaarTransform2D(
  DATATYPE* dst,            // destination
  DATATYPE* tmp,            // temporary. may be NULL
  unsigned int lenpow_x,   // log2(length_x)
  unsigned int lenpow_y)   // log2(length_y)
{
  double* t;

    if (dst==NULL)
      return ERR_GEN_NULLPOINTER;
    if ((lenpow_x<2)||(lenpow_x>16)||(lenpow_y<2)||(lenpow_y>16))
      return ERR_GEN_INVALID_PARAMETER;
    if (tmp!=NULL)
      t = tmp;
    else
      if ((t = (double*)malloc((((size_t)1)<<(lenpow_x+lenpow_y))*sizeof(double)))==NULL)
        return ERR_GEN_NOMEMORY;
    ownIPL_FastHaarTransform2D(dst,t,lenpow_x,lenpow_y);
    if (tmp==NULL)
      free(t);
    return ERR_OK;
}

#endif //DATATYPE
