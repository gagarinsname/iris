/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  20 March 2006                                          */
/*      Revision: 1.0.00                                                 */
/*      Purpose:  methods of working with ellipse                        */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/
#define __FILENUM__ 54 // __FILENUM__TAG54

#include <math.h>
#include <malloc.h>
#include "bpl.h"

// calculate inertia-equivalent ellipse parameters of a cluster of points
//   refer to BPL_ELL_GetEquivalentEllipse
int BPL_ELL_EllipseParamsOfPointGroup(
  double* ell_x,      // OUT: center position
  double* ell_y,      //
  double* ell_a,      // OUT: axis
  double* ell_b,      //
  double* ell_t,      // OUT: tilt
  const int* SelXs,   // IN:  point coordinates
  const int* SelYs,   //
  int M)              // number of points
{
  int i;
  int64 Mx,My,Mxx,Mxy,Myy;
  double tilt,coco,sisi,cosi,a,b;

    if (!M)
      return -1;
    if ((SelXs==NULL)||(SelYs==NULL))
      return -2;
    // calculate moments
    Mx = My = Mxx = Mxy = Myy = 0;
    for (i=0;i<M;i++)
    {
      Mx  += SelXs[i];
      My  += SelYs[i];
      Mxx += SelXs[i]*SelXs[i];
      Mxy += SelXs[i]*SelYs[i];
      Myy += SelYs[i]*SelYs[i];
    }
    // output first moment values
    if (ell_x)
      *ell_x = ((double)Mx)/M;
    if (ell_y)
      *ell_y = ((double)My)/M;
    // normalise second order moments
    Mxx = M*Mxx-Mx*Mx;
    Mxy = M*Mxy-Mx*My;
    Myy = M*Myy-My*My;
    if ((Mxx<=0)||(Myy<=0))
      return ERR_GEN_INVALID_PARAMETER;
    // calculate tilt, fairly simple with ATAN2
    tilt = atan2(2.*((double)Mxy),(double)(Mxx-Myy))/2.;
    if (tilt<0)
      tilt += 2*PI;
    if (tilt>PI)
      tilt -= PI;
    if (ell_t)
      *ell_t = tilt;
    //
    coco = cos(tilt);
    sisi = sin(tilt);
    cosi = 2*coco*sisi;
    coco *= coco;
    sisi *= sisi;
    a = coco*((double)Mxx)+cosi*((double)Mxy)+sisi*((double)Myy);
    b = sisi*((double)Mxx)-cosi*((double)Mxy)+coco*((double)Myy);
    if (ell_a)
      *ell_a = sqrt(a);
    if (ell_b)
      *ell_b = sqrt(b);
    return 0;
}

// interpolate group of point by circle using MSD method
// Primary formula of circle in Cartesian coordiantes is:
//   (x-x0)^2 + (y-y0)^2 = R^2
// Transforming it to form with uncorrelated linear coefficients,
// necessary for MSD, we have:
//   A*x + B*y + C = x^2 + y^2
// where A = 2*x0, B = 2*y0, C = R^2 - x0^2 - y0^2
// Thus having A,B,C found we can calculate
// x0 = A/2, y0 = B/2, R = sqrt(C + (A/2)^2 + (B/2)^2)
MIA_RESULT_CODE BPL_ELL_InterpolateByCircle(
  double* cir_x,      // OUT: center position
  double* cir_y,      //
  double* cir_r,      // OUT: axis
  double* pDiscr,     // OUT: discrepancy
  const int* SelXs,   // IN:  point coordinates
  const int* SelYs,   //
  int M)              // number of points
{
  int i,x,y,rr;
  int64 Mxx,Mxy,Myy,Mx,My,Mrrx,Mrry,Mrr;
  double adA[9],adB[3],adX[3],adS[3],XC,YC,R,D=-1.,dd;

    if ((!M)||(SelXs==NULL)||(SelYs==NULL))
      return ERR_GEN_NULLPOINTER;
    if (M<3)
      return ERR_GEN_NO_DATA;
    // calculate moments
    Mx = My = Mxx = Mxy = Myy = Mrrx = Mrry = Mrr = 0;
    for (i=0;i<M;i++)
    {
      x = SelXs[i];
      y = SelYs[i];
      rr = x*x+y*y;
      Mx   += x;
      My   += y;
      Mxx  += x*x;
      Mxy  += x*y;
      Myy  += y*y;
      Mrr  += rr;
      Mrrx += rr*x;
      Mrry += rr*y;
    }
    if ( ((double)Mxx*Myy-(double)Mxy*Mxy)/(((double)Mxx+Myy)*(Mxx+Myy)) <.00001)//   ((double)Mxx*Myy == (double)Mxy*Mxy)
      return ERR_GEN_NO_DATA;
    // fill in the data for 3rd order linear system:
    //    / Mxx Mxy Mx | Mrrx \
    //    | Mxy Myy My | Mrry |
    //    | Mx  My  M  | Mrr  /
    adA[0] = (double)Mxx;  adA[1] = (double)Mxy;  adA[2] = (double)Mx;
    adA[3] = (double)Mxy;  adA[4] = (double)Myy;  adA[5] = (double)My;
    adA[6] = (double)Mx;   adA[7] = (double)My;   adA[8] = (double)M;
    adB[0] = (double)Mrrx; adB[1] = (double)Mrry; adB[2] = (double)Mrr;
    // Solve symmetric linear system Ax=b
    BPL_Matrix_SymmetricLinearSolve(adA,adB,adX,adS,3,1);
    // calculate circle parameters
    XC = adX[0]/2.;
    YC = adX[1]/2.;
    R = sqrt(adX[2] + XC*XC + YC*YC);
    // calculate circle discrepancy (normalized on points number)
    if (pDiscr)
    {
      D = 0.;
      for (i=0;i<M;i++)
      {
        dd = sqrt((SelXs[i]-XC)*(SelXs[i]-XC)+(SelYs[i]-YC)*(SelYs[i]-YC))-R;
        if (dd<0)
          dd = -dd;
        D += dd;
      }
      D /= M;
    }
    // output data
    if (cir_x)
      *cir_x = XC;
    if (cir_y)
      *cir_y = YC;
    if (cir_r)
      *cir_r = R;
    if (pDiscr)
      *pDiscr = D;
    return 0;
}

// Interpolate group of points by ellipse using MSD method.
// Formula of ellipse in Cartesian coordinates (that we wish to get) is:
//   [(x'-x0)/a]^2 + [(y'-y0)/b]^2 = 1,                           (*)
//   where x' = x*cos(f)+y*sin(f), y' = -x*sin(f)+y*cos(f)        (**)
//   Hence, here are 5 parameters defining ellipse: x0,y0,a,b,f
// This can be rewritten as:
//   A*x^2 + B*x*y + C*y^2 + D*x + E*y = 1                        (***)
//   So A,B,C,D,E are another 5 parameters, mutually univalent to those in (*).
//   Particularly, f = atan(B/(A-C))/2
// Note that if there were no rotation, one parameters is excluded, so as to get
//   f=0 in first representation and B=0 in second one. 
// We apply MSD to formula (***) as follows: 
//   S = sum[(A*xi^2 + B*xi*yi + C*yi^2 + D*xi + E*yi - 1)^2] -> min
//   and solve linear system: { dS/dA=0, dS/dB=0, dS/dC=0, dS/dD=0, dS/dE=0 }
// Having A,B,C,D,E found, we can calculate tilt angle f. 
// Substituting it to (**): x = x'*cos(f)-y'*sin(f), y = x'*sin(f)+y'*cos(f), 
// we come to a form with excluded rotation (B=0) and uncorrelated x and y:
// A'*x'^2 + C'*y'^2 + D'*x' + E'*y' = 1
//  where A' =  A*cos^2(f) + B*cos(f)*sin(f) + C*sin^2(f)
//        C' =  A*sin^2(f) - B*cos(f)*sin(f) + C*cos^2(f)
//        D' =  D*cos(f) + E*sin(f)
//        E' = -D*sin(f) + E*cos(f)
// and from here:
//  x0 = -D'/(2*A'), y0 = -E'/(2*C')
//  a  = sqrt(A'/Q), b  = sqrt(C'/Q)
//  Q  = 1+D'^2/(4*A')+E'^2/(4*C')
MIA_RESULT_CODE BPL_ELL_InterpolateByEllipse(
  double* ell_x,      // OUT: center position
  double* ell_y,      //
  double* ell_a,      // OUT: axes
  double* ell_b,      //
  double* ell_t,      // OUT: tilt
  double* pDiscr,     // OUT: discrepancy
  const int* SelXs,   // IN:  point coordinates
  const int* SelYs,   //
  int M)              // number of points
{
  int i,x,y;
  int64 Mxxxx,Mxxxy,Mxxyy,Mxyyy,Myyyy;  // 4th order moments
  int64 Mxxx,Mxxy,Mxyy,Myyy;            // 3rd order moments
  int64 Mxx,Mxy,Myy;                    // 2nd order moments
  int64 Mx,My;                          // 1st order moments
  double adA[25],adB[5],adX[5],adS[5];
  double F,A_prim,C_prim,D_prim,E_prim,co,si,XC,YC,a,b,Q,XC_prim,YC_prim;

    if ((!M)||(SelXs==NULL)||(SelYs==NULL))
      return ERR_GEN_NULLPOINTER;
    if (M<5)
      return ERR_GEN_NO_DATA;
    // calculate moments
    Mxxxx = Mxxxy = Mxxyy = Mxyyy = Myyyy = Mx = My = 0;
    Mxxx = Mxxy = Mxyy = Myyy = Mxx = Mxy = Myy = 0;
    for (i=0;i<M;i++)
    {
      x = SelXs[i];
      y = SelYs[i];
      Mx    += x;
      My    += y;
      Mxx   += ((int64)x)*x;
      Mxy   += ((int64)x)*y;
      Myy   += ((int64)y)*y;
      Mxxx  += ((int64)x)*x*x;
      Mxxy  += ((int64)x)*x*y;
      Mxyy  += ((int64)x)*y*y;
      Myyy  += ((int64)y)*y*y;
      Mxxxx += ((int64)x)*x*x*x;
      Mxxxy += ((int64)x)*x*x*y;
      Mxxyy += ((int64)x)*x*y*y;
      Mxyyy += ((int64)x)*y*y*y;
      Myyyy += ((int64)y)*y*y*y;
    }
    if ( ((double)Mxx*Myy-(double)Mxy*Mxy)/(((double)Mxx+Myy)*(Mxx+Myy)) <.00001)
      return ERR_GEN_NO_DATA;
    // fill in the data for 5th order linear system:
    //    / Mxxxx Mxxxy Mxxyy Mxxx  Mxxy  | Mxx \
    //    | Mxxxy Mxxyy Mxyyy Mxxy  Mxyy  | Mxy |
    //    | Mxxyy Mxyyy Myyyy Mxyy  Myyy  | Myy |
    //    | Mxxx  Mxxy  Mxyy  Mxx   Mxy   | Mx  |
    //    \ Mxxy  Mxyy  Myyy  Mxy   Myy   | My  /
    adA[ 0] = (double)Mxxxx;  adA[ 1] = (double)Mxxxy;  adA[ 2] = (double)Mxxyy; adA[ 3] = (double)Mxxx; adA[ 4] = (double)Mxxy;
    adA[ 5] = (double)Mxxxy;  adA[ 6] = (double)Mxxyy;  adA[ 7] = (double)Mxyyy; adA[ 8] = (double)Mxxy; adA[ 9] = (double)Mxyy;
    adA[10] = (double)Mxxyy;  adA[11] = (double)Mxyyy;  adA[12] = (double)Myyyy; adA[13] = (double)Mxyy; adA[14] = (double)Myyy;
    adA[15] = (double)Mxxx;   adA[16] = (double)Mxxy;   adA[17] = (double)Mxyy;  adA[18] = (double)Mxx;  adA[19] = (double)Mxy;
    adA[20] = (double)Mxxy;   adA[21] = (double)Mxyy;   adA[22] = (double)Myyy;  adA[23] = (double)Mxy;  adA[24] = (double)Myy;
    adB[ 0] = (double)Mxx;    adB[ 1] = (double)Mxy;    adB[ 2] = (double)Myy;   adB[ 3] = (double)Mx;   adB[ 4] = (double)My;
    // Solve symmetric linear system Ax=b
    BPL_Matrix_SymmetricLinearSolve(adA,adB,adX,adS,5,1);
    // calculate ellipse parameters
    F = atan2(adX[1],(adX[0]-adX[2]))/2.;
    co = cos(F);
    si = sin(F);
    A_prim =  adX[0]*co*co+adX[1]*co*si+adX[2]*si*si;
    C_prim =  adX[0]*si*si-adX[1]*co*si+adX[2]*co*co;
    D_prim =  adX[3]*co+adX[4]*si;
    E_prim = -adX[3]*si+adX[4]*co;
    Q = 1+D_prim*D_prim/(4*A_prim)+E_prim*E_prim/(4*C_prim);
    // A', C' and Q should have one sign for ellipse
    if ((A_prim*Q<=0)||(C_prim*Q<=0))
      // it's not true. Some other 2nd order curve rather than ellipse
      // makes the best fit for these points.
      return ERR_GEN_INVALID_PARAMETER;
    // A', C' and Q have one sign.
    // Ellipse is the most suitable 2nd order curve.
    XC_prim = -D_prim/(2*A_prim);
    YC_prim = -E_prim/(2*C_prim);
    a = sqrt(Q/A_prim);
    b = sqrt(Q/C_prim);
    XC =  XC_prim*co-YC_prim*si;
    YC =  XC_prim*si+YC_prim*co;
    // output data
    if (ell_x)
      *ell_x = XC;
    if (ell_y)
      *ell_y = YC;
    if (ell_a)
      *ell_a = a;
    if (ell_b)
      *ell_b = b;
    if (ell_t)
      *ell_t = F;
    if (*pDiscr)
    {
      double D=0,d,xprim,yprim;
      for (i=0;i<M;i++)
      {
        xprim = (SelXs[i]-XC)*cos(F)-(SelYs[i]-YC)*sin(F);
        yprim = (SelXs[i]-XC)*sin(F)+(SelYs[i]-YC)*cos(F);
        xprim /= a;
        yprim /= b;
        d = xprim*xprim+yprim*yprim-1;
        if (d<0.)
          d = -d;
        D += d;
      }
      *pDiscr = D*sqrt(a*b)/M;
    }
    return ERR_OK;
}

// interpolate group of points by ellipse using eigen method
MIA_RESULT_CODE IPL_GEOM_InterpolateByEllipse2(
  double* ell_x,      // OUT: center position
  double* ell_y,      //
  double* ell_a,      // OUT: axes
  double* ell_b,      //
  double* ell_t,      // OUT: tilt
  double* pDiscr,     // OUT: discrepancy
  const int* SelXs,   // IN:  point coordinates
  const int* SelYs,   //
  int M)              // number of points
{
  int i,n;
  double x,y,L;
  double Mxxxx,Mxxxy,Mxxyy,Mxyyy,Myyyy;  // 4th order moments
  double Mxxx,Mxxy,Mxyy,Myyy;            // 3rd order moments
  double Mxx,Mxy,Myy;                    // 2nd order moments
  double Mx,My;                          // 1st order moments
  double adA[36],adX[36*2],adC[36],adL[6],adL2[12],adTmp[6];
  double F,A_prim,C_prim,D_prim,E_prim,co,si,XC,YC,a,b,Q,XC_prim,YC_prim;
  void* buf=NULL;
  MIA_RESULT_CODE res;

    if (!M)
      return ERR_GEN_NO_DATA;
    if ((SelXs==NULL)||(SelYs==NULL))
      return ERR_GEN_NULLPOINTER;
    // calculate moments
    Mxxxx = Mxxxy = Mxxyy = Mxyyy = Myyyy = Mx = My = 0;
    Mxxx = Mxxy = Mxyy = Myyy = Mxx = Mxy = Myy = 0;
    for (i=0;i<M;i++)
    {
      x = SelXs[i];
      y = SelYs[i];
      x /= 100;
      y /= 100;
      Mx    += x;
      My    += y;
      Mxx   += x*x;
      Mxy   += x*y;
      Myy   += y*y;
      Mxxx  += x*x*x;
      Mxxy  += x*x*y;
      Mxyy  += x*y*y;
      Myyy  += y*y*y;
      Mxxxx += x*x*x*x;
      Mxxxy += x*x*x*y;
      Mxxyy += x*x*y*y;
      Mxyyy += x*y*y*y;
      Myyyy += y*y*y*y;
    }
    // fill in 6x6 matrix:
    //    / Mxxxx Mxxxy Mxxyy Mxxx  Mxxy  Mxx \
    //    | Mxxxy Mxxyy Mxyyy Mxxy  Mxyy  Mxy |
    //    | Mxxyy Mxyyy Myyyy Mxyy  Myyy  Myy |
    //    | Mxxx  Mxxy  Mxyy  Mxx   Mxy   Mx  |
    //    | Mxxy  Mxyy  Myyy  Mxy   Myy   My  |
    //    \ Mxx   Mxy   Myy   Mx    My    M   /
    adA[ 0] = (double)Mxxxx;  adA[ 1] = (double)Mxxxy;  adA[ 2] = (double)Mxxyy; adA[ 3] = (double)Mxxx; adA[ 4] = (double)Mxxy; adA[ 5] = (double)Mxx;
    adA[ 6] = (double)Mxxxy;  adA[ 7] = (double)Mxxyy;  adA[ 8] = (double)Mxyyy; adA[ 9] = (double)Mxxy; adA[10] = (double)Mxyy; adA[11] = (double)Mxy;
    adA[12] = (double)Mxxyy;  adA[13] = (double)Mxyyy;  adA[14] = (double)Myyyy; adA[15] = (double)Mxyy; adA[16] = (double)Myyy; adA[17] = (double)Myy;
    adA[18] = (double)Mxxx;   adA[19] = (double)Mxxy;   adA[20] = (double)Mxyy;  adA[21] = (double)Mxx;  adA[22] = (double)Mxy;  adA[23] = (double)Mx;
    adA[24] = (double)Mxxy;   adA[25] = (double)Mxyy;   adA[26] = (double)Myyy;  adA[27] = (double)Mxy;  adA[28] = (double)Myy;  adA[29] = (double)My;
    adA[30] = (double)Mxx;    adA[31] = (double)Mxy;    adA[32] = (double)Myy;   adA[33] = (double)Mx;   adA[34] = (double)My;   adA[35] = (double)M;
    // fill in constraint matrix
    memset(&(adC[0]),0,sizeof(adC));
    adC[2] = 2;
    adC[7] = -1;
    adC[12] = 2;
    // solution of generalized eigen problem Ax=LBx
    buf = malloc(6*6*sizeof(double)*10);
//    res = BPL_Matrix_EigenGeneralized_Symm(
  //    adA,        // IN:  left-size matrix (A)
    //  adC,        // IN:  right-side matrix (B)
//      adX,        // OUT: eigen vectors
  //    aL,       // OUT: eigen values
    //  6,            // IN:  dimensionality
//      buf,  // IN:  temporary buffer
  //    6*6*sizeof(double)*10);        // IN:  size of temporary buffer
    n = 6;
    res = BPL_Matrix_EigenGeneralized(
      adA,
      adC,
      adX,
      adL2,
      adL,
      adTmp,
      &n,
      &i);
    free(buf);
    if (res!=ERR_OK)
      return res;
    // locate real positive finite eigen value
    n = -1;
    for (i=0;i<6;i++)
      if (adL[i]!=0)  // finite (divisor is non-zero)
        if ((adL2[i*2]!=0)&&(adL2[i*2+1]==0)) // real (Re!=0, Im=0)
          if (adL2[i*2]/adL[i]>0) // positive
          {
            if (n!=-1)
            {
              n = -2; // at least two such values - error
              break;
            }
            else
              n = i;
          }
    if (n<0)
      return n;
//n = 5;
    // eigen value / just FYI
    L = adL2[n*2]/adL[n];
    // eigen vector
    for (i=0;i<6;i++)
    {
      adX[i] = adX[n*12+i*2];
      if (adX[n*2+1]!=0)
      {
        n = -3; // eigen vector not real - error
        break;
      }
    }
    if (n<0)
      return n;
    // normalize to last coeff (F)
    if (adX[5]==0)
      return -4;  // free coeff is zero - error
    for (i=0;i<5;i++)
      adX[i] /= -adX[5];
//    // Solve symmetric linear system Ax=b
  //  BPL_Matrix_SymmetricLinearSolve(adA,adB,adX,adS,5,1);
    // calculate ellipse parameters
    F = atan2(adX[1],(adX[0]-adX[2]))/2.;
    co = cos(F);
    si = sin(F);
    A_prim =  adX[0]*co*co+adX[1]*co*si+adX[2]*si*si;
    C_prim =  adX[0]*si*si-adX[1]*co*si+adX[2]*co*co;
    D_prim =  adX[3]*co+adX[4]*si;
    E_prim = -adX[3]*si+adX[4]*co;
    Q = 1+D_prim*D_prim/(4*A_prim)+E_prim*E_prim/(4*C_prim);
    // A', C' and Q shold have one sign for ellipse
    if ((A_prim*Q<=0)||(C_prim*Q<=0))
      // it's not true. Some other 2nd order curve rather than ellipse
      // makes the best fit for these points.
      return ERR_GEN_INVALID_PARAMETER;
    // A', C' and Q have one sign.
    // Ellipse is the most suiltable 2nd order curve.
    XC_prim = -D_prim/(2*A_prim);
    YC_prim = -E_prim/(2*C_prim);
    a = sqrt(Q/A_prim);
    b = sqrt(Q/C_prim);
    XC =  XC_prim*co-YC_prim*si;
    YC =  XC_prim*si+YC_prim*co;
    // output data
    if (ell_x)
      *ell_x = XC*100;
    if (ell_y)
      *ell_y = YC*100;
    if (ell_a)
      *ell_a = a*100;
    if (ell_b)
      *ell_b = b*100;
    if (ell_t)
      *ell_t = F;
    if (*pDiscr)
      *pDiscr = 0;
    return ERR_OK;
}

// special version for ES: no doubles, no int64
MIA_RESULT_CODE BPL_ELL_InterpolateByEllipseFloat(
  float* ell_x,       // OUT: center position
  float* ell_y,       //
  float* ell_a,       // OUT: axes
  float* ell_b,       //
  float* ell_t,       // OUT: tilt
  float* pDiscr,      // OUT: discrepancy
  const int* SelXs,   // IN:  point coordinates
  const int* SelYs,   //
  int M)              // number of points
{
  int i,x,y,j;
  float Mxxxx,Mxxxy,Mxxyy,Mxyyy,Myyyy;  // 4th order moments
  float Mxxx,Mxxy,Mxyy,Myyy;            // 3rd order moments
  float Mxx,Mxy,Myy;                    // 2nd order moments
  float Mx,My;                          // 1st order moments
  int Mxxxx_p,Mxxxy_p,Mxxyy_p,Mxyyy_p,Myyyy_p;  // partial sums
  int Mxxx_p,Mxxy_p,Mxyy_p,Myyy_p;              //
  int Mxx_p,Mxy_p,Myy_p;                        //
  int Mx_p,My_p;                                //
  int maxabs,m_xx,m_xy,m_yy;
  float adA[25],adB[5],adX[5],adS[5];
  float F,A_prim,C_prim,D_prim,E_prim,co,si,XC,YC,a,b,Q,XC_prim,YC_prim;

    if ((!M)||(SelXs==NULL)||(SelYs==NULL))
      return ERR_GEN_NULLPOINTER;
    if (M<5)
      return ERR_GEN_NO_DATA;
    // estimate summation rules
    maxabs = 0;
    for (i=0;i<M;i++)
    {
      x = SelXs[i];
      if (x>=0)
      {
        if (maxabs<x)
          maxabs = x;
      }
      else
        if (maxabs<-x)
          maxabs = -x;
      y = SelYs[i];
      if (y>=0)
      {
        if (maxabs<y)
          maxabs = y;
      }
      else
        if (maxabs<-y)
          maxabs = -y;
    }
    if ((maxabs>1290)||(M>0x7fff))  // 1290=2^(31/3)
      return ERR_GEN_SIZE_NOT_MATCH;
    Mxxxx = Mxxxy = Mxxyy = Mxyyy = Myyyy = Mx = My = 0.f;
    Mxxx = Mxxy = Mxyy = Myyy = Mxx = Mxy = Myy = 0.f;
    Mxxxx_p = Mxxxy_p = Mxxyy_p = Mxyyy_p = Myyyy_p = Mx_p = My_p = 0;
    Mxxx_p = Mxxy_p = Mxyy_p = Myyy_p = Mxx_p = Mxy_p = Myy_p = 0;
    // calculate moments
    if (maxabs<=215)    // 215=2^(31/4)
    { // max value in 4th power does not exceed 0x7fffffff
      // this allows us to make all multiplications integer
      for (i=0,j=215;i<M;i++)
      {
        x = SelXs[i];
        y = SelYs[i];
        // 1st order (integer container, since 215*0x7fff<2^31)
        Mx_p    += x;
        My_p    += y;
        // 2nd order (integer container, since 215*215*0x7fff<2^31)
        Mxx_p   += (m_xx = x*x);
        Mxy_p   += (m_xy = x*y);
        Myy_p   += (m_yy = y*y);
        // 3rd order (flushing to float in 215 cycles, 215*215*215*215<2^31)
        Mxxx_p  += m_xx*x;
        Mxxy_p  += m_xx*y;
        Mxyy_p  += x*m_yy;
        Myyy_p  += y*m_yy;
        if (i==j)
        { // 215 cycles elapsed, flush 3rd order
          j = i+215;
          Mxxx += Mxxx_p;
          Mxxx_p = 0;
          Mxxy += Mxxy_p;
          Mxxy_p = 0;
          Mxyy += Mxyy_p;
          Mxyy_p = 0;
          Myyy += Myyy_p;
          Myyy_p = 0;
        }
        // 4th order ('write-thru' flushing)
        Mxxxx_p += m_xx*m_xx;
        // 0x7fffffff-215^4=10733022>57^4 - enough for buffering to make sence
        if  (Mxxxx_p>= (0x7fffffff-215*215*215*215))
        {
          Mxxxx += Mxxxx_p;
          Mxxxx_p = 0;
        }
        Mxxxy_p += m_xx*m_xy;
        if ((Mxxxy_p>= (0x7fffffff-215*215*215*215))||
            (Mxxxy_p<=-(0x7fffffff-215*215*215*215)))
        {
          Mxxxy += Mxxxy_p;
          Mxxxy_p = 0;
        }
        Mxxyy_p += m_xx*m_yy;
        if  (Mxxyy_p>= (0x7fffffff-215*215*215*215))
        {
          Mxxyy += Mxxyy_p;
          Mxxyy_p = 0;
        }
        Mxyyy_p += m_xy*m_yy;
        if ((Mxyyy_p>= (0x7fffffff-215*215*215*215))||
            (Mxyyy_p<=-(0x7fffffff-215*215*215*215)))
        {
          Mxyyy += Mxyyy_p;
          Mxyyy_p = 0;
        }
        Myyyy_p += m_yy*m_yy;
        if  (Myyyy_p>= (0x7fffffff-215*215*215*215))
        {
          Myyyy += Myyyy_p;
          Myyyy_p = 0;
        }
      }
    }
    else
    { // max value in 4th power does exceeds 0x7fffffff
      // still, max value is below 1290. 
      // Thus, all multiplication except 4th power can be integer
      for (i=0,j=1290;i<M;i++)
      {
        x = SelXs[i];
        y = SelYs[i];
        // 1st order (integer container, since 1290*0x7fff<2^31)
        Mx_p    += x;
        My_p    += y;
        // 2nd order (flushing to float in 1290 cycles, 1290*1290*1290<2^31)
        Mxx_p   += (m_xx = x*x);
        Mxy_p   += (m_xy = x*y);
        Myy_p   += (m_yy = y*y);
        if (i==j)
        { // 1290 cycles elapsed, flush 2nd order
          j = i+1290;
          Mxx += Mxx_p;
          Mxx_p = 0;
          Mxy += Mxy_p;
          Mxy_p = 0;
          Myy += Myy_p;
          Myy_p = 0;
        }
        // 3rd order ('write-thru' flushing)
        Mxxx_p  += m_xx*x;
        // 0x7fffffff-1290^3=794647>92^3 - enough for buffering to make sence
        if ((Mxxx_p>= (0x7fffffff-1290*1290*1290))||
            (Mxxx_p<=-(0x7fffffff-1290*1290*1290)))
        {
          Mxxx += Mxxx_p;
          Mxxx_p = 0;
        }
        Mxxy_p  += m_xx*y;
        if ((Mxxy_p>= (0x7fffffff-1290*1290*1290))||
            (Mxxy_p<=-(0x7fffffff-1290*1290*1290)))
        {
          Mxxy += Mxxy_p;
          Mxxy_p = 0;
        }
        Mxyy_p  += x*m_yy;
        if ((Mxyy_p>= (0x7fffffff-1290*1290*1290))||
            (Mxyy_p<=-(0x7fffffff-1290*1290*1290)))
        {
          Mxyy += Mxyy_p;
          Mxyy_p = 0;
        }
        Myyy_p  += y*m_yy;
        if ((Myyy_p>= (0x7fffffff-1290*1290*1290))||
            (Myyy_p<=-(0x7fffffff-1290*1290*1290)))
        {
          Myyy += Myyy_p;
          Myyy_p = 0;
        }
        // 4th order (no buffering, float multiplication)
        Mxxxx += ((float)m_xx)*((float)m_xx);
        Mxxxy += ((float)m_xx)*((float)m_xy);
        Mxxyy += ((float)m_xx)*((float)m_yy);
        Mxyyy += ((float)m_xy)*((float)m_yy);
        Myyyy += ((float)m_yy)*((float)m_yy);
      }
    }
    if ( ((double)Mxx_p*Myy_p-(double)Mxy_p*Mxy_p)/(((double)Mxx_p+Myy_p)*(Mxx_p+Myy_p)) <.00001)
      return ERR_GEN_NO_DATA;

    // anyway, flush all buffers
    Mx += Mx_p;
    My += My_p;
    Mxx += Mxx_p;
    Mxy += Mxy_p;
    Myy += Myy_p;
    Mxxx += Mxxx_p;
    Mxxy += Mxxy_p;
    Mxyy += Mxyy_p;
    Myyy += Myyy_p;
    Mxxxx += Mxxxx_p;
    Mxxxy += Mxxxy_p;
    Mxxyy += Mxxyy_p;
    Mxyyy += Mxyyy_p;
    Myyyy += Myyyy_p;
    // normalise moments to enhance matrix conditionality
    Mx /= maxabs;
    My /= maxabs;
    i = maxabs*maxabs;
    Mxx /= i;
    Mxy /= i;
    Myy /= i;
    F = ((float)i)*maxabs;
    Mxxx /= F;
    Mxxy /= F;
    Mxyy /= F;
    Myyy /= F;
    F = ((float)i)*i;
    Mxxxx /= F;
    Mxxxy /= F;
    Mxxyy /= F;
    Mxyyy /= F;
    Myyyy /= F;
    // fill in the data for 5th order linear system:
    //    / Mxxxx Mxxxy Mxxyy Mxxx  Mxxy  | Mxx \
    //    | Mxxxy Mxxyy Mxyyy Mxxy  Mxyy  | Mxy |
    //    | Mxxyy Mxyyy Myyyy Mxyy  Myyy  | Myy |
    //    | Mxxx  Mxxy  Mxyy  Mxx   Mxy   | Mx  |
    //    \ Mxxy  Mxyy  Myyy  Mxy   Myy   | My  /
    adA[ 0] = Mxxxx;  adA[ 1] = Mxxxy;  adA[ 2] = Mxxyy; adA[ 3] = Mxxx; adA[ 4] = Mxxy;
    adA[ 5] = Mxxxy;  adA[ 6] = Mxxyy;  adA[ 7] = Mxyyy; adA[ 8] = Mxxy; adA[ 9] = Mxyy;
    adA[10] = Mxxyy;  adA[11] = Mxyyy;  adA[12] = Myyyy; adA[13] = Mxyy; adA[14] = Myyy;
    adA[15] = Mxxx;   adA[16] = Mxxy;   adA[17] = Mxyy;  adA[18] = Mxx;  adA[19] = Mxy;
    adA[20] = Mxxy;   adA[21] = Mxyy;   adA[22] = Myyy;  adA[23] = Mxy;  adA[24] = Myy;
    adB[ 0] = Mxx;    adB[ 1] = Mxy;    adB[ 2] = Myy;   adB[ 3] = Mx;   adB[ 4] = My;
    // Solve symmetric linear system Ax=b
    BPL_Matrix_SymmetricLinearSolveFloat(adA,adB,adX,adS,5,1);
    // calculate ellipse parameters
    F  = (float)atan2(adX[1],(adX[0]-adX[2]))/2.f;
    co = (float)cos(F);
    si = (float)sin(F);
    A_prim =  adX[0]*co*co+adX[1]*co*si+adX[2]*si*si;
    C_prim =  adX[0]*si*si-adX[1]*co*si+adX[2]*co*co;
    D_prim =  adX[3]*co+adX[4]*si;
    E_prim = -adX[3]*si+adX[4]*co;
    Q = 1+D_prim*D_prim/(4*A_prim)+E_prim*E_prim/(4*C_prim);
    // A', C' and Q shold have one sign for ellipse
    if ((A_prim*Q<=0)||(C_prim*Q<=0))
      // it's not true. Some other 2nd order curve rather than ellipse
      // makes the best fit for these points.
      return ERR_GEN_INVALID_PARAMETER;
    // A', C' and Q have one sign.
    // Ellipse is the most suiltable 2nd order curve.
    XC_prim = -D_prim/(2*A_prim);
    YC_prim = -E_prim/(2*C_prim);
    a = (float)sqrt(Q/A_prim);
    b = (float)sqrt(Q/C_prim);
    XC =  XC_prim*co-YC_prim*si;
    YC =  XC_prim*si+YC_prim*co;
    // small prevention
    if (a<b)
    {
      Q = a;
      a = b;
      b = Q;
      F += PI/2;
      if (F>PI)
        F -= PI;
    }
    // output data
    if (ell_x)
      *ell_x = XC*maxabs;
    if (ell_y)
      *ell_y = YC*maxabs;
    if (ell_a)
      *ell_a = a*maxabs;
    if (ell_b)
      *ell_b = b*maxabs;
    if (ell_t)
      *ell_t = F;
    if (*pDiscr)
      *pDiscr = 0;
    return ERR_OK;
}
