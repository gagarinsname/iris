/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  13 June 2002                                           */
/*      Modified: 14 April 2003 - rotate transformation added            */
/*      Revision: 1.0.00                                                 */
/*      Purpose:  Image affine mapping by means of lest squares method   */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#define __FILENUM__ 1 // __FILENUM__TAG1

#include <memory.h>
#include <malloc.h>
#include <math.h>
#include "affine.h"
#include "ipl.h"

void SImageAffineMapping_SImageAffineMapping(
  SImageAffineMapping* pThis)
{
    pThis->S_xb_xb =
    pThis->S_xb_yb =
    pThis->S_yb_yb =
    pThis->S_xb =
    pThis->S_yb =
    pThis->S_1 =
    pThis->S_xt_xb =
    pThis->S_xt_yb =
    pThis->S_xt =
    pThis->S_yt_xb =
    pThis->S_yt_yb =
    pThis->S_yt = 0.;
    pThis->isSolved = 0;
    BPL_LIN_AllocateSystem(&(pThis->cls_x),3,3);
    BPL_LIN_AllocateSystem(&(pThis->cls_y),3,3);
}

void SImageAffineMapping_UpdateMapping(
  SImageAffineMapping* pThis,
  double xb,
  double yb,
  double xt,
  double yt)
{
    pThis->S_xb_xb += xb*xb;
    pThis->S_xb_yb += xb*yb;
    pThis->S_yb_yb += yb*yb;
    pThis->S_xb += xb;
    pThis->S_yb += yb;
    pThis->S_1 += 1.;
    pThis->S_xt_xb += xt*xb;
    pThis->S_xt_yb += xt*yb;
    pThis->S_xt += xt;
    pThis->S_yt_xb += yt*xb;
    pThis->S_yt_yb += yt*yb;
    pThis->S_yt += yt;
    pThis->isSolved = 0;
}

void SImageAffineMapping_UnupdateMapping(
  SImageAffineMapping* pThis,
  double xb,
  double yb,
  double xt,
  double yt)
{
    pThis->S_xb_xb -= xb*xb;
    pThis->S_xb_yb -= xb*yb;
    pThis->S_yb_yb -= yb*yb;
    pThis->S_xb -= xb;
    pThis->S_yb -= yb;
    pThis->S_1 -= 1.;
    pThis->S_xt_xb -= xt*xb;
    pThis->S_xt_yb -= xt*yb;
    pThis->S_xt -= xt;
    pThis->S_yt_xb -= yt*xb;
    pThis->S_yt_yb -= yt*yb;
    pThis->S_yt -= yt;
    pThis->isSolved = 0;
}

MIA_RESULT_CODE SImageAffineMapping_Map(
  SImageAffineMapping* pThis,
  double* x,
  double* y)
{
  MIA_RESULT_CODE res;
  double tmpd;

    if (pThis->S_1==0.)
      return ERR_GEN_NOT_INITIALISED;
    if (pThis->S_1<4.)
    {
      *x += (pThis->S_xt-pThis->S_xb)/pThis->S_1;
      *y += (pThis->S_yt-pThis->S_yb)/pThis->S_1;
      return ERR_OK;
    }
    if (!pThis->isSolved)
    {
      pThis->cls_x.A[0][0] = pThis->S_xb_xb;
      pThis->cls_x.A[0][1] = pThis->cls_x.A[1][0] = pThis->S_xb_yb;
      pThis->cls_x.A[1][1] = pThis->S_yb_yb;
      pThis->cls_x.A[0][2] = pThis->cls_x.A[2][0] = pThis->S_xb;
      pThis->cls_x.A[2][2] = pThis->S_1;
      pThis->cls_x.A[1][2] = pThis->cls_x.A[2][1] = pThis->S_yb;
      pThis->cls_x.b[0] = pThis->S_xt_xb;
      pThis->cls_x.b[1] = pThis->S_xt_yb;
      pThis->cls_x.b[2] = pThis->S_xt;
      pThis->cls_y.A[0][0] = pThis->S_xb_xb;
      pThis->cls_y.A[0][1] = pThis->cls_y.A[1][0] = pThis->S_xb_yb;
      pThis->cls_y.A[1][1] = pThis->S_yb_yb;
      pThis->cls_y.A[0][2] = pThis->cls_y.A[2][0] = pThis->S_xb;
      pThis->cls_y.A[2][2] = pThis->S_1;
      pThis->cls_y.A[1][2] = pThis->cls_y.A[2][1] = pThis->S_yb;
      pThis->cls_y.b[0] = pThis->S_yt_xb;
      pThis->cls_y.b[1] = pThis->S_yt_yb;
      pThis->cls_y.b[2] = pThis->S_yt;
      if ((res = BPL_LIN_SolveSystem(&(pThis->cls_x)))!=ERR_OK)
        return res;
      if ((res = BPL_LIN_SolveSystem(&(pThis->cls_y)))!=ERR_OK)
        return res;
      pThis->isSolved = 1;
    }
    tmpd = pThis->cls_x.x[0]*(*x)+pThis->cls_x.x[1]*(*y)+pThis->cls_x.x[2];
    *y =   pThis->cls_y.x[0]*(*x)+pThis->cls_y.x[1]*(*y)+pThis->cls_y.x[2];
    *x =   tmpd;
    return ERR_OK;
}

MIA_RESULT_CODE SImageAffineMapping_GetCoeffs(
  SImageAffineMapping* pThis,
  double *Coeffs)
{
    if (Coeffs==NULL)
      return ERR_GEN_NULLPOINTER;
    if (!pThis->isSolved)
      return ERR_GEN_NOT_INITIALISED;
    Coeffs[0] = pThis->cls_x.x[0];
    Coeffs[1] = pThis->cls_x.x[1];
    Coeffs[2] = pThis->cls_y.x[0];
    Coeffs[3] = pThis->cls_y.x[1];
    Coeffs[4] = pThis->cls_x.x[2];
    Coeffs[5] = pThis->cls_y.x[2];
    return ERR_OK;
}

//== rotate transformation =============================================================
typedef struct
{
  double S_xl;        // sums
  double S_xr;        //
  double S_yl;        //
  double S_yr;        //
  double S_xl2;       //
  double S_yl2;       //
  double S_xlxr;      //
  double S_xlyr;      //
  double S_ylxr;      //
  double S_ylyr;      //
  double phi;         // angle
  double k;           // scaling
  double dx;          // shift along x
  double dy;          // shift along y
  unsigned int nPt;  // number of points
  int isSolved;       // are data (phi,dx,dy) actual
} SAdjust_ImageRotateMapping, *PSAdjust_ImageRotateMapping;

// allocate structure
MIA_RESULT_CODE IPL_ADJ_CreateImageRotateMapping(
  PSAdjust_ImageRotateMapping* pps)
{
  SAdjust_ImageRotateMapping* pirm;

    if (pps==NULL)
      return ERR_GEN_NULLPOINTER;
    if ((*pps = (PSAdjust_ImageRotateMapping)(pirm = (SAdjust_ImageRotateMapping*)
            malloc(sizeof(SAdjust_ImageRotateMapping))))==NULL)
      return ERR_GEN_NOMEMORY;
    memset(pirm,0,sizeof(*pirm));
    return ERR_OK;
}

// free structure
MIA_RESULT_CODE IPL_ADJ_FreeImageRotateMapping(
  PSAdjust_ImageRotateMapping ps)
{
    if (ps==NULL)
      return ERR_GEN_NULLPOINTER;
    free(ps);
    return ERR_OK;
}

// update mapping
MIA_RESULT_CODE IPL_ADJ_UpdateImageRotateMapping(
  PSAdjust_ImageRotateMapping ps,
  double xl,
  double yl,
  double xr,
  double yr)
{
    SAdjust_ImageRotateMapping* pirm;
    if (ps==NULL)
      return ERR_GEN_NULLPOINTER;
    pirm = (SAdjust_ImageRotateMapping*)ps;
    pirm->S_xl += xl;
    pirm->S_yl += yl;
    pirm->S_xr += xr;
    pirm->S_yr += yr;
    pirm->S_xl2 += xl*xl;       //
    pirm->S_yl2 += yl*yl;       //
    pirm->S_xlxr += xl*xr;      //
    pirm->S_xlyr += xl*yr;      //
    pirm->S_ylxr += yl*xr;      //
    pirm->S_ylyr += yl*yr;      //
    pirm->nPt++;
    pirm->isSolved = 0;
    return ERR_OK;
}

// un-update mapping
MIA_RESULT_CODE IPL_ADJ_UnupdateImageRotateMapping(
  PSAdjust_ImageRotateMapping ps,
  double xl,
  double yl,
  double xr,
  double yr)
{
    SAdjust_ImageRotateMapping* pirm = (SAdjust_ImageRotateMapping*)ps;
    if (ps==NULL)
      return ERR_GEN_NULLPOINTER;
    if (!pirm->nPt)
      return ERR_GEN_NOT_INITIALISED;
    pirm->S_xl -= xl;
    pirm->S_yl -= yl;
    pirm->S_xr -= xr;
    pirm->S_yr -= yr;
    pirm->S_xl2 -= xl*xl;       //
    pirm->S_yl2 -= yl*yl;       //
    pirm->S_xlxr -= xl*xr;      //
    pirm->S_xlyr -= xl*yr;      //
    pirm->S_ylxr -= yl*xr;      //
    pirm->S_ylyr -= yl*yr;      //
    pirm->isSolved = 0;
    pirm->nPt--;
    return ERR_OK;
}

// calculate rotation mapping
staticFunc MIA_RESULT_CODE ownIPL_ADJ_SolveImageRotateMapping(
  SAdjust_ImageRotateMapping* pirm)
{
    if (pirm->nPt==0)
      return ERR_GEN_NOT_INITIALISED;
    if (pirm->nPt>=4)
    { /*
      pirm->phi = atan2((pirm->S_xl*pirm->S_yr - pirm->nPt*pirm->S_xlyr) -
                        (pirm->S_xr*pirm->S_yl - pirm->nPt*pirm->S_ylxr),
                        (pirm->S_xl*pirm->S_xr - pirm->nPt*pirm->S_xlxr) +
                        (pirm->S_yl*pirm->S_yr - pirm->nPt*pirm->S_ylyr));
      pirm->dx = (pirm->S_xr-cos(pirm->phi)*pirm->S_xl+sin(pirm->phi)*pirm->S_yl)/pirm->nPt;
      pirm->dy = (pirm->S_yr-sin(pirm->phi)*pirm->S_xl-cos(pirm->phi)*pirm->S_yl)/pirm->nPt;
      */
    double phi = pirm->phi = atan2((pirm->nPt*pirm->S_xlyr - pirm->S_xl*pirm->S_yr) -
                        (pirm->nPt*pirm->S_ylxr - pirm->S_xr*pirm->S_yl),
                        (pirm->nPt*pirm->S_xlxr - pirm->S_xl*pirm->S_xr) +
                        (pirm->nPt*pirm->S_ylyr - pirm->S_yl*pirm->S_yr));
      pirm->k = (cos(phi)*pirm->nPt*pirm->S_xlxr - sin(phi)*pirm->nPt*pirm->S_xlyr -
                 cos(phi)*pirm->S_xl*pirm->S_xr + sin(phi)*pirm->S_xl*pirm->S_yr +
                 cos(phi)*pirm->nPt*pirm->S_ylyr + sin(phi)*pirm->nPt*pirm->S_ylxr -
                 cos(phi)*pirm->S_yl*pirm->S_yr - sin(phi)*pirm->S_yl*pirm->S_xr) /
                (pirm->nPt*pirm->S_xl2 - pirm->S_xl*pirm->S_xl +
                 pirm->nPt*pirm->S_yl2 - pirm->S_yl*pirm->S_yl);
      pirm->dx = (pirm->S_xr-pirm->k*cos(pirm->phi)*pirm->S_xl+pirm->k*sin(pirm->phi)*pirm->S_yl)/pirm->nPt;
      pirm->dy = (pirm->S_yr-pirm->k*sin(pirm->phi)*pirm->S_xl-pirm->k*cos(pirm->phi)*pirm->S_yl)/pirm->nPt;
    }
    else
    {
      pirm->k = 1.;
      pirm->phi = 0.;
      pirm->dx = (pirm->S_xr-pirm->S_xl)/pirm->nPt;
      pirm->dx = (pirm->S_yr-pirm->S_yl)/pirm->nPt;
    }
    pirm->isSolved = 1;
    return ERR_OK;
}

// get mapping coefficients
MIA_RESULT_CODE IPL_ADJ_GetImageRotateMapping(
  PSAdjust_ImageRotateMapping ps,
  double* pphi,
  double* pdx,
  double* pdy)
{
  SAdjust_ImageRotateMapping* pirm;
  MIA_RESULT_CODE res;

    if (ps==NULL)
      return ERR_GEN_NULLPOINTER;
    if (!((pirm = (SAdjust_ImageRotateMapping*)ps)->isSolved))
      if ((res = ownIPL_ADJ_SolveImageRotateMapping(pirm))!=ERR_OK)
        return res;
    *pphi = pirm->phi;
    *pdx = pirm->dx;
    *pdy = pirm->dy;
    return ERR_OK;
}
