/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  19 March 2002                                          */
/*      Revision: 1.0.00                                                 */
/*      Purpose:  Fast Walsh transform                                   */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#ifdef DATATYPE
#include "ipl.h"

// perform one butterfly level of Walsh transform
staticFunc void ownIPL_FWT_Stage(
  DATATYPE* p,
  unsigned int aperture,
  unsigned int size)
{
    while (size--)
      if (size&aperture)
      {
        p[size-aperture] = p[size-aperture]+  p[size];  // a=a0+b0
        p[size] =          p[size-aperture]-2*p[size];  // b=a0-b0=(a0+b0)-2*b0=a-2*b0
      }
}

// transforming 1-dimensional data with power of 2 length
// no checks are made since function should work as fast as possible
// however this is not actual since used only inside 2D tranform, 
// which does all the checks
staticFunc void ownIPL_FastWalshTransform1D(
  DATATYPE* dst,          // destination
  unsigned int size)     // length
{
  DATATYPE t;
  unsigned int i,j,m;

    // butterflies
    for (i=size>>1;i;i>>=1)
      ownIPL_FWT_Stage(dst,i,size);
    // bit reverse (adapted from Numerical Recipies fft)
    // Anchestor algorithm was in Fortran thus starting indexing from 1. 
    // Adjust to C's zero-indexing. Doing one this small operation is really 
    // better way than sweatful rewriting of this toughly designed code. 
    dst--;
    for (i=1,j=1;i<size;i++)
    {
      if (j>i)
      {
        t = dst[i];
        dst[i] = dst[j];
        dst[j] = t;
      }
      m = size>>1;
      while (m &&(j>m))
      {
        j -= m;
        m >>= 1;
      }
      j += m;
    }
}

// transforming 2-dimensional data with power of 2 sizes
staticFunc void ownIPL_FastWalshTransform2D(
  DATATYPE* dst,            // destination
  unsigned int size_x,     // length along x
  unsigned int size_y)     // length along y
{
  unsigned int i,j,m,k;
  DATATYPE t;

    // butterflies for y and for x
    for (i=size_x*size_y>>1;i;i>>=1)
      ownIPL_FWT_Stage(dst,i,size_x*size_y);
    // bit reverse for x
    dst--;
    for (i=1,j=1;i<size_x;i++)
    {
      if (j>i)
        for (k=0;k<size_y;k++)
        {
          t = dst[i+k*size_x];
          dst[i+k*size_x] = dst[j+k*size_x];
          dst[j+k*size_x] = t;
        }
      m = size_x>>1;
      while (m &&(j>m))
      {
        j -= m;
        m >>= 1;
      }
      j += m;
    }
    // bit reverse for y
    dst -= (size_x-1);
    for (i=1,j=1;i<size_y;i++)
    {
      if (j>i)
        for (k=0;k<size_x;k++)
        {
          t = dst[i*size_x+k];
          dst[i*size_x+k] = dst[j*size_x+k];
          dst[j*size_x+k] = t;
        }
      m = size_y>>1;
      while (m &&(j>m))
      {
        j -= m;
        m >>= 1;
      }
      j += m;
    }
}

// 1D Walsh transform of DATATYPE
MIA_RESULT_CODE IPL_FastWalshTransform1D(
  DATATYPE* dst,          // destination
  unsigned int lenpow)   // log2(length)
{
    if (dst==NULL)
      return ERR_GEN_NULLPOINTER;
    if ((lenpow<2)||(lenpow>16))
      return ERR_GEN_INVALID_PARAMETER;
    ownIPL_FastWalshTransform1D(dst,1<<lenpow);
    return ERR_OK;
}

// 2D Walsh transform of DATATYPE
MIA_RESULT_CODE IPL_FastWalshTransform2D(
  DATATYPE* dst,            // destination
  unsigned int lenpow_x,   // log2(length_x)
  unsigned int lenpow_y)   // log2(length_y)
{
    if (dst==NULL)
      return ERR_GEN_NULLPOINTER;
    if ((lenpow_x<2)||(lenpow_x>16)||(lenpow_y<2)||(lenpow_y>16))
      return ERR_GEN_INVALID_PARAMETER;
    ownIPL_FastWalshTransform2D(dst,1<<lenpow_x,1<<lenpow_y);
    return ERR_OK;
}

#endif //DATATYPE
