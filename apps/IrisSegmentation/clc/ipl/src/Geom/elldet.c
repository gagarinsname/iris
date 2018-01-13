/*------------------------------------------------------------------------*/
/*                                                                        */
/*      Created:  20 September 2014                                       */
/*      Revision: 1.0.00                                                  */
/*      Purpose:  ellipse detection                                       */
/*      Authors:                                                          */
/*        Ivan Matveev                                                    */
/*                                                                        */
/*------------------------------------------------------------------------*/

#define __FILENUM__ 142 // __FILENUM__TAG142

#include <stdio.h>
#include <malloc.h>
#include <memory.h>
#include "bpl.h"
#include "ipl.h"

// get the best ellipse
MIA_RESULT_CODE IPL_ELL_Detect_v1(
  double* pdElls,   // OUT: 5 ellipse parameters - x,y,a,b,f
  const uint8* im,  // IN:  source image
  int xs,           // IN:  image size
  int ys)
{
  int idx,nObjnum;
  uint8* binim,binlev;
  SConObjPrimaryParams* pCOPP;
  double moments[6],ellpar[5],bestfillqual=0.,fillqual;
  MIA_RESULT_CODE res=ERR_OK;
  
    // clear
    memset(pdElls,0,5*sizeof(*pdElls));
    // alloc
    binim = (uint8*)malloc(xs*ys);
    pCOPP = (SConObjPrimaryParams*)malloc(xs*ys);
    for (binlev=10;binlev<=30;binlev+=5)
    {
      // binarize
      for (idx=xs*ys-1;idx>=0;idx--)
        binim[idx] = (im[idx]>binlev)?0:255;
  //DBGL_FILE_SaveUint8Image(binim,xs,xs,ys,"c:\\1bibi.bmp");
      // conobj
      nObjnum = xs*ys/sizeof(*pCOPP);
      if ((res = IPL_CONOB_EnumerateConnectedObjects(
          pCOPP,  // OUT: buffer for object parameters
          &nObjnum,         // IN/OUT: number of objects allocated/found
          binim,    // input image (completely zeroed inside the function)
          xs,              // image stride in elements
          xs,               // image size
          ys,               //
          50,      // lower object mass threshold (-1 - disable)
          300,     // higher object mass threshold (-1 - disable)
          0))!=ERR_OK)        // if(TRUE) connect 8 neibors else - 4
        goto Detect_v1_exit;
      for (idx=0;idx<nObjnum;idx++)
      {
        // copy moments
        moments[0] = pCOPP[idx].M;
        moments[1] = pCOPP[idx].MX;
        moments[2] = pCOPP[idx].MY;
        moments[3] = (double)(pCOPP[idx].MXX);
        moments[4] = (double)(pCOPP[idx].MXY);
        moments[5] = (double)(pCOPP[idx].MYY);
        // calculate inertia-equivalent ellipse parameters from simple moments
        res = BPL_ELL_GetEquivalentEllipse(
          ellpar,         // OUT: xc,yc,a,b,tilt
          moments); // IN:  M,Mx,My,Mxx,Mxy,Myy
        if (res!=ERR_OK)
          continue;
        if (ellpar[2]>2.*ellpar[3])
          continue;
        fillqual = moments[0]/(PI*ellpar[2]*ellpar[3]);
        if ((fillqual>bestfillqual)&&(fillqual>.7))
        {
          bestfillqual = fillqual;
          memcpy(pdElls,ellpar,5*sizeof(*pdElls));
        }
      }
    }
Detect_v1_exit:;
    if (binim)
      free(binim);
    if (pCOPP)
      free(pCOPP);
    return res;
}

// detect multiple ellipses
MIA_RESULT_CODE IPL_ELL_DetectMulti_v1(
  SEllipseExt** ppsElls,  // OUT: set of detected ellipses with qualities
  int* pnElls,            // OUT: number of detected ellipses
  const uint8* im,        // IN:  source image
  int xs,                 // IN:  image size
  int ys)
{
  int idx,nObjnum;
  uint8* binim;
  SConObjPrimaryParams* pCOPP;
  double moments[6],ellpar[5],fillqual;
  MIA_RESULT_CODE res=ERR_OK;
  SEllipseExt* psElls=NULL;
  int nell_a=0,nell_u=0;
  int cThr,validx,aTvals[] = {10, 15, 21, 28, 36, 45};
  
    // alloc
    binim = (uint8*)malloc(xs*ys);
    pCOPP = (SConObjPrimaryParams*)malloc(xs*ys);
    for (validx=0;validx<sizeof(aTvals)/sizeof(*aTvals);validx++)
    {
      // binarize
      cThr = aTvals[validx];
      for (idx=xs*ys-1;idx>=0;idx--)
        binim[idx] = (im[idx]>cThr)?0:255;
  //DBGL_FILE_SaveUint8Image(binim,xs,xs,ys,"c:\\1bibi.bmp");
      // conobj
      nObjnum = xs*ys/sizeof(*pCOPP);
      if ((res = IPL_CONOB_EnumerateConnectedObjects(
          pCOPP,  // OUT: buffer for object parameters
          &nObjnum,         // IN/OUT: number of objects allocated/found
          binim,    // input image (completely zeroed inside the function)
          xs,              // image stride in elements
          xs,               // image size
          ys,               //
          50,      // lower object mass threshold (-1 - disable)
          300,     // higher object mass threshold (-1 - disable)
          0))!=ERR_OK)        // if(TRUE) connect 8 neibors else - 4
        goto DetectMulti_v1_exit;
      for (idx=0;idx<nObjnum;idx++)
      {
        // copy moments
        moments[0] = pCOPP[idx].M;
        moments[1] = pCOPP[idx].MX;
        moments[2] = pCOPP[idx].MY;
        moments[3] = (double)(pCOPP[idx].MXX);
        moments[4] = (double)(pCOPP[idx].MXY);
        moments[5] = (double)(pCOPP[idx].MYY);
        // calculate inertia-equivalent ellipse parameters from simple moments
        res = BPL_ELL_GetEquivalentEllipse(
          ellpar,         // OUT: xc,yc,a,b,tilt
          moments); // IN:  M,Mx,My,Mxx,Mxy,Myy
        if (res!=ERR_OK)
          continue;
        if (ellpar[2]>2.*ellpar[3])
          continue;
        fillqual = moments[0]/(PI*ellpar[2]*ellpar[3]);
        if (fillqual<.7)
          continue;
        // memorize the ellipse
        if (nell_a==nell_u)
          psElls = (SEllipseExt*)realloc(psElls,(nell_a+=1024)*sizeof(*psElls));
        memcpy(psElls[nell_u].ellpar,ellpar,5*sizeof(*ellpar));
        psElls[nell_u].Qfill = fillqual;
        psElls[nell_u].Qelong = ellpar[3]/ellpar[2];
        psElls[nell_u].Qdir = (1.+ellpar[4]*ellpar[4])/2.;
        psElls[nell_u].binval = cThr;
        psElls[nell_u].multi = 1;
        nell_u++;
      }
    }
DetectMulti_v1_exit:;
    if (binim)
      free(binim);
    if (pCOPP)
      free(pCOPP);
    *ppsElls = psElls;
    *pnElls = nell_u;
    return res;
}

// detect multiple ellipses
MIA_RESULT_CODE IPL_ELL_DetectMultiExt_v1(
  SEllipseExt** ppsElls,  // OUT: set of detected ellipses with qualities
  int* pnElls,            // OUT: number of detected ellipses
  const uint8* im,        // IN:  source image
  int xs,                 // IN:  image size
  int ys,
  const SEllDetPar* pEp)  // IN:  parameters of ellipse selection
{
  int idx,nObjnum;
  uint8* binim;
  SConObjPrimaryParams* pCOPP;
  double moments[6],ellpar[5],fillqual;
  MIA_RESULT_CODE res=ERR_OK;
  SEllipseExt* psElls=NULL;
  int nell_a=0,nell_u=0;
  int cThr,validx;
  
    // alloc
    binim = (uint8*)malloc(xs*ys);
    pCOPP = (SConObjPrimaryParams*)malloc(xs*ys);
    for (validx=0;validx<pEp->nThresholds;validx++)
    {
      // binarize
      cThr = pEp->aThresholds[validx];
      if (cThr>0)
        for (idx=xs*ys-1;idx>=0;idx--)
          binim[idx] = (im[idx]>cThr)?0:255;
      else
        for (idx=xs*ys-1;idx>=0;idx--)
          binim[idx] = (im[idx]>-cThr)?255:0;
//DBGL_FILE_SaveUint8Image(binim,xs,xs,ys,"c:\\1bibi.bmp");
      // conobj
      nObjnum = xs*ys/sizeof(*pCOPP);
      if ((res = IPL_CONOB_EnumerateConnectedObjects(
          pCOPP,  // OUT: buffer for object parameters
          &nObjnum,         // IN/OUT: number of objects allocated/found
          binim,    // input image (completely zeroed inside the function)
          xs,              // image stride in elements
          xs,               // image size
          ys,               //
          pEp->nMassMin,      // lower object mass threshold (-1 - disable)
          pEp->nMassMax,     // higher object mass threshold (-1 - disable)
          0))!=ERR_OK)        // if(TRUE) connect 8 neibors else - 4
        goto DetectMulti_v2_exit;
      for (idx=0;idx<nObjnum;idx++)
      {
        // copy moments
        moments[0] = pCOPP[idx].M;
        moments[1] = pCOPP[idx].MX;
        moments[2] = pCOPP[idx].MY;
        moments[3] = (double)(pCOPP[idx].MXX);
        moments[4] = (double)(pCOPP[idx].MXY);
        moments[5] = (double)(pCOPP[idx].MYY);
        // calculate inertia-equivalent ellipse parameters from simple moments
        res = BPL_ELL_GetEquivalentEllipse(
          ellpar,         // OUT: xc,yc,a,b,tilt
          moments); // IN:  M,Mx,My,Mxx,Mxy,Myy
        if (res!=ERR_OK)
          continue;
        if (ellpar[3]/ellpar[2]<pEp->QelongMin)
          continue;
        fillqual = moments[0]/(PI*ellpar[2]*ellpar[3]);
        if (fillqual<pEp->QfillMin)
          continue;
        // memorize the ellipse
        if (nell_a==nell_u)
          psElls = (SEllipseExt*)realloc(psElls,(nell_a+=1024)*sizeof(*psElls));
        memcpy(psElls[nell_u].ellpar,ellpar,5*sizeof(*ellpar));
        psElls[nell_u].Qfill = fillqual;
        psElls[nell_u].Qelong = ellpar[3]/ellpar[2];
        psElls[nell_u].Qdir = (1.+ellpar[4]*ellpar[4])/2.;
        psElls[nell_u].binval = cThr;
        psElls[nell_u].multi = 1;
        nell_u++;
      }
    }
DetectMulti_v2_exit:;
    if (binim)
      free(binim);
    if (pCOPP)
      free(pCOPP);
    *ppsElls = psElls;
    *pnElls = nell_u;
    return res;
}
