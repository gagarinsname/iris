/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:     4 July 2001                                         */
/*      Revision:    1.0.00                                              */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#include <math.h>
#include <malloc.h>
#include "../ownbpl.h"

MIA_RESULT_CODE BPL_LIN_AllocateSystem(
  SLinearSystem* pThis,
  int nStr,
  int nCol)
{
  int i;

    if ((pThis->A = (double**)malloc(nStr*sizeof(double*)+
                                     nStr*nCol*sizeof(double)+
                                     2*nCol*sizeof(double)))==NULL)
      return ERR_GEN_NOMEMORY;
    pThis->A[0] = (double*)(pThis->A+nStr);
    for (i=1;i<nStr;i++)
      pThis->A[i] = pThis->A[i-1]+nCol;
    pThis->b = pThis->A[nStr-1]+nCol;
    pThis->x = pThis->b+nStr;
    pThis->isSolved = 0;
    pThis->nCol = nCol;
    pThis->nStr = nStr;
    return ERR_OK;
}

void BPL_LIN_FreeSystem(
  SLinearSystem* pThis)
{
    if (pThis->A)
      free(pThis->A);
    pThis->A = NULL;
}

MIA_RESULT_CODE BPL_LIN_SolveSystem(
  SLinearSystem* pThis)
{
  unsigned int i,j,l;
  double SqLength,LengthA,LengthQ,Scalar;
  unsigned int Size;
  double** A1 = NULL;   // Matrix A1 =  U(l)*..*U(1)* {MatrA|VectB}   !! Size x (Size+1)
  double* omega = NULL;
  MIA_RESULT_CODE ret = ERR_OK;

    if (pThis->A==NULL)
      return ERR_GEN_NOT_INITIALISED;
    if (pThis->nCol!=pThis->nStr)
      return MIA_LIN_MATRIX_NOT_SQUARE;
    Size = pThis->nCol;
    for (;;)
    {
      // Allocating memory
        ret = ERR_GEN_NOMEMORY;
      if ( ( omega = (double*)malloc(Size*sizeof(double))) ==NULL)
        break;
      if ( ( A1 = (double**)malloc(Size*sizeof(double))) ==NULL)
        break;
      if ( ( A1[0] = (double*)malloc(Size*(Size+1)*sizeof(double))) ==NULL)
        break;
      for (i=1; i<Size; i++)
        A1[i] = A1[i-1] + Size+1;
      for (i=0; i<Size; i++)                      // A1 = {MatrA|VectB}
        for (j=0; j<Size; j++)
          A1[i][j] = pThis->A[i][j];
      for (i=0; i<Size; i++)
          A1[i][Size] = pThis->b[i];
        ret = ERR_OK;
      for (l=0; l<Size-1; l++)        // Loop of iteration
      {
        i=l+1;                                   //
        while ( (i<Size) ? (A1[i][l]==0) : 0 )   // if the rest of l'th column is
          i++;                                   // parallel with (1,0,...0) then
        if (i<Size)                              //  do nothing
        {
          SqLength = 0.;
          for (i=l+1; i<Size; i++)            // Calculating of omega :
            SqLength += A1[i][l]*A1[i][l];    // \/
          LengthA = sqrt ( SqLength + A1[l][l]*A1[l][l] );
          if (A1[l][l]>0)
          {
            omega[l] = A1[l][l] + LengthA;                        // (1)
            A1[l][l] = -LengthA;
            LengthQ = sqrt ( SqLength + omega[l]*omega[l] );
            omega[l] /= LengthQ;
          }
          else
          {
            omega[l] = A1[l][l] - LengthA;
            A1[l][l] = +LengthA;                                  // (2)
            LengthQ = sqrt ( SqLength + omega[l]*omega[l] );
            omega[l] /= LengthQ;
          }
          for (i=l+1; i<Size; i++)                          //
            omega[i] = A1[i][l] / LengthQ;                  //
          //              | E   0  |
          //      newA1 = | 0  U(l)| * A1 ; U(l) = E - 2*omega*transp(omega)
          //
          for (i=l+1; i<Size; i++)        // y = e*Length(y);  y  - rest of l'th column
            A1[i][l] = 0.;                // First component was calculated in (1) or (2)
          for (j=l+1; j<=Size; j++)       // newh = h - Scalar*omega ;
          {                               //         Scalar = 2*(omega,h)
            Scalar = omega[l]*A1[l][j];
            for (i=l+1; i<Size; i++)
              Scalar += omega[i]*A1[i][j];
            Scalar *= 2;
            for (i=l; i<Size; i++)
              A1[i][j] -= Scalar * omega[i];
          }
        }
      }
      //    Now we have A1 in the form :    / * * * . * | * \
      //                                    | 0 * * . * | * |
      //                                    | 0 0 * . * | * |
      //                                    | . . . . * |   |
      //                                    \ 0 0 0 0 * | * /
      //  We can easy solve it from end
      //
      for (l=Size-1; l<Size; l--)
      {
        Scalar = A1[l][Size];
        for (i=l+1; i<Size; i++)
          Scalar -= A1[l][i] * pThis->x[i];
        if (A1[l][l]==0)
        {
          ret = MIA_LIN_MATRIX_DEGENERATE;
          break;
        }
        pThis->x[l] = Scalar / A1[l][l];
      }
      break;
    }
    // Freeing Memory
    if (omega)
      free (omega);
    if (A1)
    {
      if (A1[0])
        free(A1[0]);
      free(A1);
    }
    return ret;
}

// calculate inertia-equivalent ellipse parameters from simple moments
MIA_RESULT_CODE BPL_ELL_GetEquivalentEllipse(
  double* ellpar,         // OUT: xc,yc,a,b,tilt
  const double* moments)  // IN:  M,Mx,My,Mxx,Mxy,Myy
{
  double tilt,coco,sisi,cosi,a,b,Mxx,Mxy,Myy;

    if ((moments==NULL)||(ellpar==NULL))
      return ERR_GEN_NULLPOINTER;
    if (moments[0]<=0)
      return ERR_GEN_INVALID_PARAMETER;
    // output first moment values
    ellpar[0] = moments[1]/moments[0];
    ellpar[1] = moments[2]/moments[0];
    // normalize second order moments
    Mxx = (moments[3]-moments[1]*moments[1]/moments[0])/moments[0];
    Mxy = (moments[4]-moments[1]*moments[2]/moments[0])/moments[0];
    Myy = (moments[5]-moments[2]*moments[2]/moments[0])/moments[0];
    // calculate tilt, fairly simple with ATAN2
    tilt = atan2(2.*Mxy,Mxx-Myy)/2.;
    if (tilt<0)
      tilt += 2*PI;
    if (tilt>PI)
      tilt -= PI;
    ellpar[4] = tilt;
    //
    coco = cos(tilt);
    sisi = sin(tilt);
    cosi = 2*coco*sisi;
    coco *= coco;
    sisi *= sisi;
    a = coco*((double)Mxx)+cosi*((double)Mxy)+sisi*((double)Myy);
    b = sisi*((double)Mxx)-cosi*((double)Mxy)+coco*((double)Myy);
    ellpar[2] = 2*sqrt(a);
    ellpar[3] = 2*sqrt(b);
    // small prevention
    if (ellpar[2]<ellpar[3])
    {
      coco = ellpar[2];
      ellpar[2] = ellpar[3];
      ellpar[3] = coco;
      ellpar[4] += PI/2;
      if (ellpar[4]>PI)
        ellpar[4] -= PI;
    }
/*    // normalize
    a = sqrt(moments[0]/(PI*ellpar[2]*ellpar[3]));
    ellpar[2] *= a;
    ellpar[3] *= a;*/
    return ERR_OK;
}
