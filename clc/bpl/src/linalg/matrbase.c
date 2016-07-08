/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:     13 April 2001                                       */
/*      Modified:    14 December 2004                                    */
/*      Revision:    2.0.00                                              */
/*      Purpose:     Basic matrix operations                             */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#include <math.h>
#include "bpl.h"

// transpose
void BPL_Matrix_Transpose(
  double* matr,
  int size)
{
  int i,j;
  double tmpvar;

    for (i=0;i<size;i++)
      for (j=0;j<i;j++)
      {
        tmpvar = matr[i*size+j];
        matr[i*size+j] = matr[j*size+i];
        matr[j*size+i] = tmpvar;
      }
}

MIA_RESULT_CODE BPL_Matrix_Multiply(
  double *R,        // OUT: result matrix       dimensions depend on mode
  const double *A,  // IN:  first multiplier    height:width = l:m
  const double *B,  // IN:  second multiplier   height:width = m:n
  int aw,           // IN:  matrix A width
  int ah,           // IN:  matrix A height
  int bw,           // IN:  matrix B width
  int bh,           // IN:  matrix B height
  int mode)         // IN:  multiplication mode
{
  int i,j,k;
  double s;

    // check arguments
    if ((R==NULL)||(A==NULL)||(B==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((aw<=0)||(ah<=0)||(bw<=0)||(bh<=0))
      return ERR_GEN_NO_DATA;
    // four possible variants
    switch (mode)
    {
      case IPL_MATRIX_MULTIPLY_AB:
        if (aw!=bh)
          return ERR_GEN_SIZE_NOT_MATCH;
        for (i=0;i<ah;i++)
          for (j=0;j<bw;j++)
          {
            s = 0.;
            for (k=0;k<aw;k++)
              s += A[i*aw+k]*B[k*bw+j];
            R[i*bw+j] = s;
          }
        return ERR_OK;
      case IPL_MATRIX_MULTIPLY_ABt:
        if (aw!=bw)
          return ERR_GEN_SIZE_NOT_MATCH;
        for (i=0;i<ah;i++)
          for (j=0;j<bh;j++)
          {
            s = 0.;
            for (k=0;k<aw;k++)
              s += A[i*aw+k]*B[j*bw+k];
            R[i*bh+j] = s;
          }
        return ERR_OK;
      case IPL_MATRIX_MULTIPLY_AtB:
        if (ah!=bh)
          return ERR_GEN_SIZE_NOT_MATCH;
        for (i=0;i<aw;i++)
          for (j=0;j<bw;j++)
          {
            s = 0.;
            for (k=0;k<ah;k++)
              s += A[k*aw+i]*B[k*bw+j];
            R[i*bw+j] = s;
          }
        return ERR_OK;
      case IPL_MATRIX_MULTIPLY_AtBt:
        if (ah!=bw)
          return ERR_GEN_SIZE_NOT_MATCH;
        for (i=0;i<aw;i++)
          for (j=0;j<bh;j++)
          {
            s = 0.;
            for (k=0;k<ah;k++)
              s += A[k*aw+i]*B[j*bw+k];
            R[i*bh+j] = s;
          }
        return ERR_OK;
    }
    return ERR_GEN_INVALID_PARAMETER;
}

