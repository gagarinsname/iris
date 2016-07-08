/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:     14 December 2004                                    */
/*      Revision:    1.0.00                                              */
/*      Purpose:     Checking matrix                                     */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#include <math.h>
#include <malloc.h>
#include "bpl.h"

MIA_RESULT_CODE BPL_Matrix_EstimateErrorNonUnitar(
  const double* A,
  int n,
  double* pdErrMaxAbsDevDiag,   // OUT: maximum absolute deviation of diagonal elements from 1
  double* pdErrMaxAbsDevNdiag,  // OUT: maximum absolute deviation of off-diagonal elements from 0
  double* pdErrAvgAbsDevDiag,   // OUT: average absolute deviation of diagonal elements from 1
  double* pdErrAvgAbsDevNdiag)  // OUT: average absolute deviation of off-diagonal elements from 0
{
  int i,j;
  double ErrMaxAbsDevDiag;   // maximum absolute deviation of diagonal elements from 1
  double ErrMaxAbsDevNdiag;  // maximum absolute deviation of off-diagonal elements from 0
  double ErrAvgAbsDevDiag;   // average absolute deviation of diagonal elements from 1
  double ErrAvgAbsDevNdiag;  // average absolute deviation of off-diagonal elements from 0
  double* pB;

    // check args
    if (A==NULL)
      return ERR_GEN_NULLPOINTER;
    if (n<=0)
      return ERR_GEN_INVALID_PARAMETER;
    // allocate
    if ((pB = (double*)malloc(n*n*sizeof(pB[0])))==NULL)
      return ERR_GEN_NOMEMORY;
    // multiple
    BPL_Matrix_Multiply(pB,A,A,n,n,n,n,IPL_MATRIX_MULTIPLY_ABt);
    // check how close is to zero
    ErrMaxAbsDevDiag = ErrMaxAbsDevNdiag =
    ErrAvgAbsDevDiag = ErrAvgAbsDevNdiag = 0.;
    for (i=0;i<n;i++)
    {
      if (ErrMaxAbsDevDiag<fabs(pB[i*n+i]-1.))
        ErrMaxAbsDevDiag = fabs(pB[i*n+i]-1.);
      ErrAvgAbsDevDiag += fabs(pB[i*n+i]-1.);
      for (j=0;j<n;j++)
        if (i!=j)
        {
          if (ErrMaxAbsDevNdiag<fabs(pB[i*n+j]))
            ErrMaxAbsDevNdiag = fabs(pB[i*n+j]);
          ErrAvgAbsDevNdiag += fabs(pB[i*n+j]);
        }
    }
    // output
    if (pdErrMaxAbsDevDiag)
      *pdErrMaxAbsDevDiag = ErrMaxAbsDevDiag ;
    if (pdErrMaxAbsDevNdiag)
      *pdErrMaxAbsDevNdiag = ErrMaxAbsDevNdiag;
    if (pdErrAvgAbsDevDiag)
      *pdErrAvgAbsDevDiag = ErrAvgAbsDevDiag/n;
    if (pdErrAvgAbsDevNdiag)
      *pdErrAvgAbsDevNdiag = ErrAvgAbsDevNdiag/(n*(n-1));
    // Free resources
    free(pB);
    return ERR_OK;
}

MIA_RESULT_CODE BPL_Matrix_EstimateErrorNonEigen(
  const double* A,    // IN:  source matrix
  const double* Vec,  // IN:  eigen vectors
  const double* Val,  // IN:  eigen values
  int n,              // IN:  size of source matrix
  int m,              // IN:  number of eigen vectors to be tested
  double* maxdiff)
{
  int i,j,k;
  double s,dev,maxdev;

    // check args
    if ((A==NULL)||(Vec==NULL)||(Val==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((n<=0)||(m<=0))
      return ERR_GEN_INVALID_PARAMETER;
    //
    maxdev = 0.;
    for (i=0;i<m;i++)       // vectors
    {
      dev = 0.;
      for (j=0;j<n;j++)     // i-th vector components
      {
        // calculate j-th component of i-th vector
        s = 0.;
        for (k=0;k<n;k++)
          s += A[j*n+k]*Vec[i*n+k];
        // now s = A*x(j), will be s = A*x(j)-l(i)*x(j)
        s -= Val[i]*Vec[i*n+j];
        // add to total deviation
        dev += s*s;
      }
      dev = sqrt(dev/Val[i]);
      if (dev>maxdev)
        maxdev = dev;
    }
    if (maxdiff)
      *maxdiff = maxdev;
    return ERR_OK;
}

// test AQAt=E
MIA_RESULT_CODE BPL_Matrix_EstimateErrorNonUnitarGeneral(
  const double* A,
  const double* Q,
  int n,
  double* pdErrMaxAbsDevDiag,   // OUT: maximum absolute deviation of diagonal elements from 1
  double* pdErrMaxAbsDevNdiag,  // OUT: maximum absolute deviation of off-diagonal elements from 0
  double* pdErrAvgAbsDevDiag,   // OUT: average absolute deviation of diagonal elements from 1
  double* pdErrAvgAbsDevNdiag)  // OUT: average absolute deviation of off-diagonal elements from 0
{
  int i,j;
  double ErrMaxAbsDevDiag;   // maximum absolute deviation of diagonal elements from 1
  double ErrMaxAbsDevNdiag;  // maximum absolute deviation of off-diagonal elements from 0
  double ErrAvgAbsDevDiag;   // average absolute deviation of diagonal elements from 1
  double ErrAvgAbsDevNdiag;  // average absolute deviation of off-diagonal elements from 0
  double *pB,*pC;

    // check args
    if ((A==NULL)||(Q==NULL))
      return ERR_GEN_NULLPOINTER;
    if (n<=0)
      return ERR_GEN_INVALID_PARAMETER;
    // allocate
    if ((pB = (double*)malloc(2*n*n*sizeof(pB[0])))==0)
      return -1;
    pC = pB+n*n;
    // multiple
    BPL_Matrix_Multiply(pC,A,Q,n,n,n,n,IPL_MATRIX_MULTIPLY_AB);
    BPL_Matrix_Multiply(pB,pC,A,n,n,n,n,IPL_MATRIX_MULTIPLY_ABt);
    // check how close is to zero
    ErrMaxAbsDevDiag = ErrMaxAbsDevNdiag =
    ErrAvgAbsDevDiag = ErrAvgAbsDevNdiag = 0.;
    for (i=0;i<n;i++)
    {
      if (ErrMaxAbsDevDiag<fabs(pB[i*n+i]-1.))
        ErrMaxAbsDevDiag = fabs(pB[i*n+i]-1.);
      ErrAvgAbsDevDiag += fabs(pB[i*n+i]-1.);
      for (j=0;j<n;j++)
        if (i!=j)
        {
          if (ErrMaxAbsDevNdiag<fabs(pB[i*n+j]))
            ErrMaxAbsDevNdiag = fabs(pB[i*n+j]);
          ErrAvgAbsDevNdiag += fabs(pB[i*n+j]);
        }
    }
    // output
    if (pdErrMaxAbsDevDiag)
      *pdErrMaxAbsDevDiag = ErrMaxAbsDevDiag ;
    if (pdErrMaxAbsDevNdiag)
      *pdErrMaxAbsDevNdiag = ErrMaxAbsDevNdiag;
    if (pdErrAvgAbsDevDiag)
      *pdErrAvgAbsDevDiag = ErrAvgAbsDevDiag/n;
    if (pdErrAvgAbsDevNdiag)
      *pdErrAvgAbsDevNdiag = ErrAvgAbsDevNdiag/(n*(n-1));
    // Free resources
    free(pB);
    return 0;
}

// Ax=LBx
MIA_RESULT_CODE BPL_Matrix_EstimateErrorNonEigenGeneral(
  const double* A,    // IN:  source matrix at left
  const double* B,    // IN:  source matrix at right
  const double* Vec,  // IN:  eigen vectors
  const double* Val,  // IN:  eigen values
  int n,              // IN:  size of source matrix
  int m,              // IN:  number of eigen vectors to be tested
  double* maxdiff)
{
  int i,j,k;
  double sl,dev,maxdev,sr;

    // check args
    if ((A==NULL)||(B==NULL)||(Vec==NULL)||(Val==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((n<=0)||(m<=0))
      return ERR_GEN_INVALID_PARAMETER;
    maxdev = 0.;
    for (i=0;i<m;i++)       // vectors
    {
      dev = 0.;
      for (j=0;j<n;j++)     // i-th vector components
      {
        // calculate j-th component of i-th vector
        sl = sr = 0.;
        for (k=0;k<n;k++)
        {
          sl += A[j*n+k]*Vec[i*n+k];
          sr += B[j*n+k]*Vec[i*n+k];
        }
        // now sl = A*x(j), sr = B*x(j)
        sl -= Val[i]*sr;
        // add to total deviation
        dev += sl*sl;
      }
      dev = sqrt(dev/Val[i]);
      if (dev>maxdev)
        maxdev = dev;
    }
    if (maxdiff)
      *maxdiff = maxdev;
    return 0;
}
