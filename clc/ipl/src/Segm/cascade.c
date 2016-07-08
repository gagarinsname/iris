/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  31 December 2013                                       */
/*      Revision: 1.0.00                                                 */
/*      Purpose:  Viola-Jones own implementation (quite poor)            */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/
#define __FILENUM__ 156 // __FILENUM__TAG156

#include <malloc.h>
#include <stdio.h>
#include <memory.h>
#include "errorcodes.h"
#include "stddefs.h"
#include "ipl.h"
//#include "dbg.h"
#include "../../../../api/inc/cgl.h"

#define THM 1

// primitive types are given by bit mask
//  0/1 - horizontal/vertical;
//  0/2 - 2-sided/3-sided
//  0/4 - positive/negative

#define VJ_PRI_HORIZONTAL 0
#define VJ_PRI_VERTICAL   1
#define VJ_PRI_2SIDED     0
#define VJ_PRI_3SIDED     2
#define VJ_PRI_POSITIVE   0
#define VJ_PRI_NEGATIVE   4

static MIA_RESULT_CODE ownVJ_CreatePrimitives(
  SViolaJonesPrimitive** ppPri,
  int* pnPriCount,
  int xs,
  int ys,
  int str)
{
  MIA_RESULT_CODE res=ERR_OK;
  SViolaJonesPrimitive* pPri = NULL;
  int x,y,w,h,nPri_u=0,nPri_a=0;
/*
    if (nPri_u==nPri_a)
      if ((pPri = (SViolaJonesPrimitive*)realloc(pPri,sizeof(*pPri)*(nPri_a+=1024)))==NULL)
      {
        res = ERR_GEN_NOMEMORY;
        goto ownVJ_CreatePrimitives_exit;
      }
    memset(pPri+nPri_u,0,sizeof(*pPri));
    // horizontal positive
    pPri[nPri_u].nType = VJ_PRI_HORIZONTAL|VJ_PRI_2SIDED;
    pPri[nPri_u].idx[0] = 0;
    pPri[nPri_u].idx[1] = xs/2;
    pPri[nPri_u].idx[2] = xs;
    pPri[nPri_u].idx[3] = ys*str;
    pPri[nPri_u].idx[4] = ys*str+xs/2;
    pPri[nPri_u].idx[5] = ys*str+xs;
    nPri_u++;
goto ownVJ_CreatePrimitives_exit;
*/
    // -- 2-sided horizontal
    for (w=2;w<=xs;w+=2)
      for (h=1;h<=ys;h++)
        for (x=0;x<=xs-w;x++)
          for (y=0;y<=ys-h;y++)
          {
            if (nPri_u==nPri_a)
              if ((pPri = (SViolaJonesPrimitive*)realloc(pPri,sizeof(*pPri)*(nPri_a+=1024)))==NULL)
              {
                res = ERR_GEN_NOMEMORY;
                goto ownVJ_CreatePrimitives_exit;
              }
            memset(pPri+nPri_u,0,sizeof(*pPri));
            // horizontal positive
            pPri[nPri_u].nType = VJ_PRI_HORIZONTAL|VJ_PRI_2SIDED;
            pPri[nPri_u].idx[0] =  y   *str+x;
            pPri[nPri_u].idx[1] =  y   *str+x+w/2;
            pPri[nPri_u].idx[2] =  y   *str+x+w;
            pPri[nPri_u].idx[3] = (y+h)*str+x;
            pPri[nPri_u].idx[4] = (y+h)*str+x+w/2;
            pPri[nPri_u].idx[5] = (y+h)*str+x+w;
            nPri_u++;
          }
    // -- 2-sided vertical
    for (w=1;w<=xs;w++)
      for (h=2;h<=ys;h+=2)
        for (x=0;x<=xs-w;x++)
          for (y=0;y<=ys-h;y++)
          {
            if (nPri_u==nPri_a)
              if ((pPri = (SViolaJonesPrimitive*)realloc(pPri,sizeof(*pPri)*(nPri_a+=1024)))==NULL)
              {
                res = ERR_GEN_NOMEMORY;
                goto ownVJ_CreatePrimitives_exit;
              }
            memset(pPri+nPri_u,0,sizeof(*pPri));
            // vertical positive
            pPri[nPri_u].nType = VJ_PRI_VERTICAL|VJ_PRI_2SIDED;
            pPri[nPri_u].idx[0] =  y     *str+x;
            pPri[nPri_u].idx[1] = (y+h/2)*str+x;
            pPri[nPri_u].idx[2] = (y+h  )*str+x;
            pPri[nPri_u].idx[3] =  y     *str+x+w;
            pPri[nPri_u].idx[4] = (y+h/2)*str+x+w;
            pPri[nPri_u].idx[5] = (y+h  )*str+x+w;
            nPri_u++;
          }
    // -- 3-sided horizontal
    for (w=4;w<=xs;w+=4)
      for (h=1;h<=ys;h++)
        for (x=0;x<=xs-w;x++)
          for (y=0;y<=ys-h;y++)
          {
            if (nPri_u==nPri_a)
              if ((pPri = (SViolaJonesPrimitive*)realloc(pPri,sizeof(*pPri)*(nPri_a+=1024)))==NULL)
              {
                res = ERR_GEN_NOMEMORY;
                goto ownVJ_CreatePrimitives_exit;
              }
            memset(pPri+nPri_u,0,sizeof(*pPri));
            // horizontal positive
            pPri[nPri_u].nType = VJ_PRI_HORIZONTAL|VJ_PRI_3SIDED;
            pPri[nPri_u].idx[0] =  y   *str+x;
            pPri[nPri_u].idx[1] =  y   *str+x+w/4;
            pPri[nPri_u].idx[2] =  y   *str+x+3*w/4;
            pPri[nPri_u].idx[3] =  y   *str+x+w;
            pPri[nPri_u].idx[4] = (y+h)*str+x;
            pPri[nPri_u].idx[5] = (y+h)*str+x+w/4;
            pPri[nPri_u].idx[6] = (y+h)*str+x+3*w/4;
            pPri[nPri_u].idx[7] = (y+h)*str+x+w;
            nPri_u++;
          }
    // -- 3-sided vertical
    for (w=1;w<=xs;w++)
      for (h=4;h<=ys;h+=4)
        for (x=0;x<=xs-w;x++)
          for (y=0;y<=ys-h;y++)
          {
            if (nPri_u==nPri_a)
              if ((pPri = (SViolaJonesPrimitive*)realloc(pPri,sizeof(*pPri)*(nPri_a+=1024)))==NULL)
              {
                res = ERR_GEN_NOMEMORY;
                goto ownVJ_CreatePrimitives_exit;
              }
            memset(pPri+nPri_u,0,sizeof(*pPri));
            // vertical positive
            pPri[nPri_u].nType = VJ_PRI_VERTICAL|VJ_PRI_3SIDED;
            pPri[nPri_u].idx[0] =  y       *str+x;
            pPri[nPri_u].idx[1] = (y+h/4)  *str+x;
            pPri[nPri_u].idx[2] = (y+3*h/4)*str+x;
            pPri[nPri_u].idx[3] = (y+h)    *str+x;
            pPri[nPri_u].idx[4] =  y       *str+x+w;
            pPri[nPri_u].idx[5] = (y+h/4)  *str+x+w;
            pPri[nPri_u].idx[6] = (y+3*h/4)*str+x+w;
            pPri[nPri_u].idx[7] = (y+h)    *str+x+w;
            nPri_u++;
          }
ownVJ_CreatePrimitives_exit:;
    if (res!=ERR_OK)
    {
      if (pPri)
        free(pPri);
      pPri = NULL;
      nPri_u = 0;
    }
    *ppPri = pPri;
    *pnPriCount = nPri_u;
    return res;
}

#define DRAWBAR(C) \
for (y=y_;y<y_+h_;y++) \
  for (x=x_;x<x_+w_;x++) \
    im[y*xs+x] = C;
/*
void ownVJ_DrawPrimitive(
  const char* filename,
  SViolaJonesPrimitive* pP,
  int xs,
  int ys)
{
  uint8* im,v1,v2;
  int x_,y_,w_,h_,x,y;

    im = (uint8*)malloc(xs*ys);
    memset(im,127,xs*ys);
    if (pP->nType&VJ_PRI_NEGATIVE)
    {
      v1 = 254;
      v2 = 0;
    }
    else
    {
      v1 = 0;
      v2 = 254;
    }
    switch (pP->nType&(VJ_PRI_3SIDED|VJ_PRI_VERTICAL))
    {
      case VJ_PRI_2SIDED|VJ_PRI_HORIZONTAL:
        x_ = pP->idx[0]%(xs+1);
        y_ = pP->idx[0]/(xs+1);
        w_ = pP->idx[1]-pP->idx[0];
        h_ = (pP->idx[3]-pP->idx[0])/(xs+1);
        DRAWBAR(v1)
        x_ = pP->idx[1]%(xs+1);
        y_ = pP->idx[1]/(xs+1);
        w_ = pP->idx[2]-pP->idx[1];
        DRAWBAR(v2)
        break;
      case VJ_PRI_2SIDED|VJ_PRI_VERTICAL:
        x_ = pP->idx[0]%(xs+1);
        y_ = pP->idx[0]/(xs+1);
        w_ = pP->idx[3]-pP->idx[0];
        h_ = (pP->idx[1]-pP->idx[0])/(xs+1);
        DRAWBAR(v1)
        x_ = pP->idx[1]%(xs+1);
        y_ = pP->idx[1]/(xs+1);
        h_ = (pP->idx[2]-pP->idx[1])/(xs+1);
        DRAWBAR(v2)
        break;
      case VJ_PRI_3SIDED|VJ_PRI_HORIZONTAL:
        x_ = pP->idx[0]%(xs+1);
        y_ = pP->idx[0]/(xs+1);
        w_ = pP->idx[1]-pP->idx[0];
        h_ = (pP->idx[4]-pP->idx[0])/(xs+1);
        DRAWBAR(v1)
        x_ = pP->idx[1]%(xs+1);
        y_ = pP->idx[1]/(xs+1);
        w_ = pP->idx[2]-pP->idx[1];
        DRAWBAR(v2)
        x_ = pP->idx[2]%(xs+1);
        y_ = pP->idx[2]/(xs+1);
        w_ = pP->idx[3]-pP->idx[2];
        DRAWBAR(v1)
        break;
      case VJ_PRI_3SIDED|VJ_PRI_VERTICAL:
        x_ = pP->idx[0]%(xs+1);
        y_ = pP->idx[0]/(xs+1);
        w_ = pP->idx[4]-pP->idx[0];
        h_ = (pP->idx[1]-pP->idx[0])/(xs+1);
        DRAWBAR(v1)
        x_ = pP->idx[1]%(xs+1);
        y_ = pP->idx[1]/(xs+1);
        h_ = (pP->idx[2]-pP->idx[1])/(xs+1);
        DRAWBAR(v2)
        x_ = pP->idx[2]%(xs+1);
        y_ = pP->idx[2]/(xs+1);
        h_ = (pP->idx[3]-pP->idx[2])/(xs+1);
        DRAWBAR(v1)
        break;
    }
    DBGL_FILE_SaveUint8Image(im,xs,xs,ys,filename);
    free(im);
}
*/
MIA_RESULT_CODE VJ_train(
  SViolaJonesClassifier** ppSVJ,
  const uint8* pPosSamples,
  int nPosSamples,
  const uint8* pNegSamples,
  int nNegSamples,
  int xs,
  int ys,
  const char* repa,
  const char* nambeg)
{
  SViolaJonesClassifier* pVJ=NULL;
  MIA_RESULT_CODE res = ERR_OK;
  int cSample;
  int *pIImPos=NULL,*pIImNeg=NULL;
  SViolaJonesPrimitive* pPri=NULL;
  int nPri_u=0;
  int str,stry;
  int cCascade,cPrimitive,cStage;
  int nNegSamples_current,nBestPrim;
  double *pWeights=NULL,dBestVal,dValPos,dValNeg;
  SViolaJonesPrimitive* pri;
  SViolaJonesCascade* cas;
  int *iim,tmpint,bBestIsNeg=0;
  int nCascades_a=0;
  int nPrimitives_a=0,nNegSamples_o,nThrPos,nThrNeg;
//  char nama[FILENAME_MAX];

    // check arguments
    if ((ppSVJ==NULL)||(pPosSamples==NULL)||(pNegSamples==NULL))
      return ERR_GEN_NULLPOINTER;
    pVJ = *ppSVJ;
    if (pVJ!=NULL)
      return ERR_GEN_INITIALISED_ALREADY;
    if ((xs<=0)||(ys<=0))
      return ERR_GEN_EMPTY_ROI;
    if ((nPosSamples<=0)||(nNegSamples<=0))
      return ERR_GEN_NO_DATA;
    // allocate structure
    if ((pVJ = (SViolaJonesClassifier*)malloc(sizeof(*pVJ)))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto VJ_train_exit;
    }
    memset(pVJ,0,sizeof(*pVJ));
    // fill in
    pVJ->xs = xs;
    pVJ->ys = ys;
    str = xs+1;
    stry = ys+1;
    // create primitives
    if ((res = ownVJ_CreatePrimitives(&pPri,&nPri_u,xs,ys,str))!=ERR_OK)
      goto VJ_train_exit;
    if (repa)
      FILE_Printf(repa,"%d primitives created\n",nPri_u);
    // allocate integral images
    if ((pIImPos = (int*)malloc(sizeof(*pIImPos)*str*stry*(nPosSamples+nNegSamples)))==NULL)
    {
      if (repa)
        FILE_Printf(repa,"cannot allocate integral image");
      res = ERR_GEN_NOMEMORY;
      goto VJ_train_exit;
    }
    pIImNeg = pIImPos+str*stry*nPosSamples;
    // create integral images
    for (cSample=0;cSample<nPosSamples;cSample++)
      IPL_PMP_MakeIntegralImage_prim(
        pIImPos+str*stry*cSample,  // OUT: integer image
        pPosSamples+xs*ys*cSample,  // IN:  source image
        xs,ys);                   // IN:  size
    for (cSample=0;cSample<nNegSamples;cSample++)
      IPL_PMP_MakeIntegralImage_prim(
        pIImNeg+str*stry*cSample,  // OUT: integer image
        pNegSamples+xs*ys*cSample,  // IN:  source image
        xs,ys);                   // IN:  size
    if (repa)
      FILE_Printf(repa,"integral images filled in\n");
    if ((pWeights = (double*)malloc(sizeof(*pWeights)*(nPosSamples+nNegSamples)))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto VJ_train_exit;
    }
    // main cycle - by cascades
    nNegSamples_current = nNegSamples;
    for (cCascade=0;cCascade<7;cCascade++)
    {
      if (nNegSamples_current<=0)
        break;  // all negatives dropped, no more classification possible
      //== new cascade begun, initialize it
      if (cCascade==nCascades_a)
        if ((pVJ->pCascades = realloc(pVJ->pCascades,sizeof(pVJ->pCascades[0])*(nCascades_a+=1)))==NULL)
        {
          res = ERR_GEN_NOMEMORY;
          goto VJ_train_exit;
        }
      cas = pVJ->pCascades+cCascade;
      memset(cas,0,sizeof(*cas));
      // init weights
      for (cSample=0;cSample<nPosSamples;cSample++)
        pWeights[cSample] = .5/nPosSamples;
      for (cSample=0;cSample<nNegSamples_current;cSample++)
        pWeights[cSample+nPosSamples] = .5/nNegSamples_current;
      // cycle by stages
      nPrimitives_a = 0;
      for (cStage=0;cStage<cCascade+2;cStage++)
      {
        if (cStage==nPrimitives_a)
          if ((cas->pPrimitives = (SViolaJonesPrimitive*)realloc(cas->pPrimitives,sizeof(cas->pPrimitives[0])*(nPrimitives_a+=1)))==NULL)
          {
            res = ERR_GEN_NOMEMORY;
            goto VJ_train_exit;
          }
        // apply classifiers (i.e. primitives) to samples
        nBestPrim = 0;
        dBestVal = 1e100;
        // cycle by primitives
        for (cPrimitive=0;cPrimitive<nPri_u;cPrimitive++)
        {
          pri = pPri+cPrimitive;
          // init thresholds for positive and negative classifiers
          nThrPos =  0x7fffffff;
          nThrNeg = -0x7fffffff;
          // cycle by positive samples
          for (cSample=0;cSample<nPosSamples;cSample++)
          {
            iim = pIImPos+str*stry*cSample;
            if ((pri->nType&(VJ_PRI_3SIDED|VJ_PRI_2SIDED))==VJ_PRI_3SIDED)
              tmpint = +iim[pri->idx[0]]-2*iim[pri->idx[1]]+2*iim[pri->idx[2]]-iim[pri->idx[3]]
                       -iim[pri->idx[4]]+2*iim[pri->idx[5]]-2*iim[pri->idx[6]]+iim[pri->idx[7]];
            else
              tmpint = +iim[pri->idx[0]]-2*iim[pri->idx[1]]+iim[pri->idx[2]]
                       -iim[pri->idx[3]]+2*iim[pri->idx[4]]-iim[pri->idx[5]];
            if (nThrPos>=tmpint)
              nThrPos = tmpint-1;
            if (nThrNeg<=tmpint)
              nThrNeg = tmpint+1;
          }
          pri->nThreshold  = nThrPos;
          pri->nThreshold2 = nThrNeg;
          // init total error of positive and negative classifiers
          dValPos = dValNeg = 0.;
          // cycle by negative samples
          for (cSample=0;cSample<nNegSamples_current;cSample++)
          {
            iim = pIImNeg+str*stry*cSample;
            if ((pri->nType&(VJ_PRI_3SIDED|VJ_PRI_2SIDED))==VJ_PRI_3SIDED)
              tmpint = +iim[pri->idx[0]]-2*iim[pri->idx[1]]+2*iim[pri->idx[2]]-iim[pri->idx[3]]
                       -iim[pri->idx[4]]+2*iim[pri->idx[5]]-2*iim[pri->idx[6]]+iim[pri->idx[7]];
            else
              tmpint = +iim[pri->idx[0]]-2*iim[pri->idx[1]]+iim[pri->idx[2]]
                       -iim[pri->idx[3]]+2*iim[pri->idx[4]]-iim[pri->idx[5]];
            if (tmpint-pri->nThreshold>0)
              // positive response, accumulate error for positive
              dValPos += pWeights[cSample+nPosSamples];
            if (tmpint-pri->nThreshold2<0)
              // negative response, accumulate error for negative
              dValNeg += pWeights[cSample+nPosSamples];
          }
          // for given primitive all samples passed, select best by error
          if (dBestVal>dValPos)
          {
            dBestVal = dValPos;
            nBestPrim = cPrimitive;
            bBestIsNeg = 0;
          }
          if (dBestVal>dValNeg)
          {
            dBestVal = dValNeg;
            nBestPrim = cPrimitive;
            bBestIsNeg = 1;
          }
        }
        // passed all primitives
        if (repa)
          FILE_Printf(repa,"Cascade %d, Stage %d, Error = %f on primitive %d%c\n",
            cCascade,cStage,(float)dBestVal,nBestPrim,(bBestIsNeg?'-':'+'));
        // remember the weak classifier (i.e. primitive)
        memcpy(cas->pPrimitives+cStage,pPri+nBestPrim,sizeof(*pPri));
        if (bBestIsNeg)
        {
          cas->pPrimitives[cStage].nType |= VJ_PRI_NEGATIVE;
          cas->pPrimitives[cStage].nThreshold = cas->pPrimitives[cStage].nThreshold2;
        }
        if (dBestVal>0)
          cas->pPrimitives[cStage].weight = log((1.-dBestVal)/dBestVal);
        else
          cas->pPrimitives[cStage].weight = 100.;
#if 0
if (nambeg)
{
  sprintf(nama,"%s_%d_%d_%d.bmp",nambeg,cCascade,cStage,nBestPrim);
  ownVJ_DrawPrimitive(nama,pPri+nBestPrim,xs,ys);
}
#endif
        // update weights / for positive samples
        pri = cas->pPrimitives+cStage;
        for (cSample=0;cSample<nPosSamples;cSample++)
        {
          iim = pIImPos+str*stry*cSample;
          if ((pri->nType&(VJ_PRI_3SIDED|VJ_PRI_2SIDED))==VJ_PRI_3SIDED)
            tmpint = +iim[pri->idx[0]]-2*iim[pri->idx[1]]+2*iim[pri->idx[2]]-iim[pri->idx[3]]
                     -iim[pri->idx[4]]+2*iim[pri->idx[5]]-2*iim[pri->idx[6]]+iim[pri->idx[7]];
          else
            tmpint = +iim[pri->idx[0]]-2*iim[pri->idx[1]]+iim[pri->idx[2]]
                     -iim[pri->idx[3]]+2*iim[pri->idx[4]]-iim[pri->idx[5]];
          tmpint -= THM*pri->nThreshold;
          if ( ((tmpint>0)&&(!(pri->nType&VJ_PRI_NEGATIVE))) ||
               ((tmpint<0)&&( (pri->nType&VJ_PRI_NEGATIVE))) ) 
            // correct! classification
            pWeights[cSample] *= (dBestVal/(1.-dBestVal));
        }
        // update weights / for negative samples
        for (cSample=0;cSample<nNegSamples_current;cSample++)
        {
          iim = pIImNeg+str*stry*cSample;
          if ((pri->nType&(VJ_PRI_3SIDED|VJ_PRI_2SIDED))==VJ_PRI_3SIDED)
            tmpint = +iim[pri->idx[0]]-2*iim[pri->idx[1]]+2*iim[pri->idx[2]]-iim[pri->idx[3]]
                     -iim[pri->idx[4]]+2*iim[pri->idx[5]]-2*iim[pri->idx[6]]+iim[pri->idx[7]];
          else
            tmpint = +iim[pri->idx[0]]-2*iim[pri->idx[1]]+iim[pri->idx[2]]
                     -iim[pri->idx[3]]+2*iim[pri->idx[4]]-iim[pri->idx[5]];
          tmpint -= THM*pri->nThreshold;
          if ( ((tmpint<=0)&&(!(pri->nType&VJ_PRI_NEGATIVE))) ||
               ((tmpint>=0)&&( (pri->nType&VJ_PRI_NEGATIVE))) ) 
            // correct! classification
            pWeights[cSample+nPosSamples] *= (dBestVal/(1.-dBestVal));
        }
        // normalize weights
        dBestVal = 0.;
        for (cSample=0;cSample<nPosSamples;cSample++)
          dBestVal += pWeights[cSample];
        for (cSample=0;cSample<nNegSamples_current;cSample++)
          dBestVal += pWeights[cSample+nPosSamples];
        for (cSample=0;cSample<nPosSamples;cSample++)
          pWeights[cSample] /= dBestVal;
        for (cSample=0;cSample<nNegSamples_current;cSample++)
          pWeights[cSample+nPosSamples] /= dBestVal;
      }
      // stages ended. Weak classifiers (primitives) for Adaboost selected. 
      // record number of primitives
      cas->nPrimitives = cStage;
      // see how many positive samples are rejected (optional)
      nNegSamples_o = 0;
      for (cSample=0;cSample<nPosSamples;cSample++)
      {
        dBestVal = dValPos = 0.;
        iim = pIImPos+str*stry*cSample;
        //++ this is a strong classifier
        for (cPrimitive=0;cPrimitive<cas->nPrimitives;cPrimitive++)
        {
          pri = cas->pPrimitives+cPrimitive;
          if ((pri->nType&(VJ_PRI_3SIDED|VJ_PRI_2SIDED))==VJ_PRI_3SIDED)
            tmpint = +iim[pri->idx[0]]-2*iim[pri->idx[1]]+2*iim[pri->idx[2]]-iim[pri->idx[3]]
                     -iim[pri->idx[4]]+2*iim[pri->idx[5]]-2*iim[pri->idx[6]]+iim[pri->idx[7]];
          else
            tmpint = +iim[pri->idx[0]]-2*iim[pri->idx[1]]+iim[pri->idx[2]]
                     -iim[pri->idx[3]]+2*iim[pri->idx[4]]-iim[pri->idx[5]];
          tmpint -= THM*pri->nThreshold;
          if ( ((tmpint>0)&&(!(pri->nType&VJ_PRI_NEGATIVE))) ||
               ((tmpint<0)&&( (pri->nType&VJ_PRI_NEGATIVE))) ) 
            dBestVal += pri->weight;
          dValPos += pri->weight;
        }
        //-- strong classifier
        if (dBestVal*2<=dValPos)
          nNegSamples_o++;
      }
      if (repa)
        FILE_Printf(repa,"%d positive samples are erroneously rejected\n",nNegSamples_o);
      // run strong classifier and throw away rejected negative samples
      nNegSamples_o = nNegSamples_current;
      for (cSample=0;cSample<nNegSamples_current;cSample++)
      {
        dBestVal = dValPos = 0.;
        iim = pIImNeg+str*stry*cSample;
        //++ this is a strong classifier
        for (cPrimitive=0;cPrimitive<cas->nPrimitives;cPrimitive++)
        {
          pri = cas->pPrimitives+cPrimitive;
          if ((pri->nType&(VJ_PRI_3SIDED|VJ_PRI_2SIDED))==VJ_PRI_3SIDED)
            tmpint = +iim[pri->idx[0]]-2*iim[pri->idx[1]]+2*iim[pri->idx[2]]-iim[pri->idx[3]]
                     -iim[pri->idx[4]]+2*iim[pri->idx[5]]-2*iim[pri->idx[6]]+iim[pri->idx[7]];
          else
            tmpint = +iim[pri->idx[0]]-2*iim[pri->idx[1]]+iim[pri->idx[2]]
                     -iim[pri->idx[3]]+2*iim[pri->idx[4]]-iim[pri->idx[5]];
          tmpint -= THM*pri->nThreshold;
          if ( ((tmpint>0)&&(!(pri->nType&VJ_PRI_NEGATIVE))) ||
               ((tmpint<0)&&( (pri->nType&VJ_PRI_NEGATIVE))) ) 
            dBestVal += pri->weight;
          dValPos += pri->weight;
        }
        //-- strong classifier
        if (dBestVal*2<=dValPos)
        { // throw image away
          memmove(iim,pIImNeg+str*stry*(nNegSamples_current-1),str*stry*sizeof(*iim));
          nNegSamples_current--;
          cSample--;
        }
      }
      if (repa)
        FILE_Printf(repa,"%d negative samples dropped, %d left\n",nNegSamples_o-nNegSamples_current,nNegSamples_current);
    }
    pVJ->nCascade = cCascade;
VJ_train_exit:;
    if (pPri)
      free(pPri);
    if (pWeights)
      free(pWeights);
    if (pIImPos)
      free(pIImPos);
    if (res!=ERR_OK)
      VJ_Free(&pVJ);
    *ppSVJ = pVJ;
    return res;
}

MIA_RESULT_CODE VJ_Save(
  const SViolaJonesClassifier* pS,
  const char* nama)
{
  int cCascade,cPrimitive;
  SViolaJonesCascade* pC;
  SViolaJonesPrimitive* pP;

    if ((pS==NULL)||(nama==NULL))
      return ERR_GEN_NULLPOINTER;
    if (pS->pCascades==NULL)
      return ERR_GEN_NULLPOINTER;
    FILE_Printf(nama,"%d %d %d\n",pS->xs,pS->ys,pS->nCascade);
    for (cCascade=0;cCascade<pS->nCascade;cCascade++)
    {
      pC = pS->pCascades+cCascade;
      FILE_Printf(nama,"  %d\n",pC->nPrimitives);
//      if ((pC->pWeights==NULL)||(pC->pPrimitives==NULL))
  //      return ERR_GEN_NULLPOINTER;
      if (pC->pPrimitives==NULL)
        return ERR_GEN_NULLPOINTER;
      for (cPrimitive=0;cPrimitive<pC->nPrimitives;cPrimitive++)
      {
        pP = pC->pPrimitives+cPrimitive;
        FILE_Printf(nama,"    %f    %d %d    %d %d %d %d %d %d %d %d\n",
          (float)(pP->weight),    pP->nType,pP->nThreshold,
          pP->idx[0],pP->idx[1],pP->idx[2],pP->idx[3],
          pP->idx[4],pP->idx[5],pP->idx[6],pP->idx[7]);
//        FILE_Printf(nama,"    %f    %d %d    %d %d %d %d    %d %d %d %d\n",
  //        (float)(pC->pWeights),    pP->nType,pP->nThreshold,
    //      pP->x[0],pP->x[1],pP->x[2],pP->x[3],pP->y[0],pP->y[1],pP->y[2],pP->y[3]);
      }
    }
    return ERR_OK;
}

MIA_RESULT_CODE VJ_Load(
  const SViolaJonesClassifier** ppS,
  const char* nama)
{
  int len,idx=0,xs,ys;
  char* buf=NULL;
  int cCascade,cPrimitive;
  SViolaJonesClassifier* pS=NULL;
  SViolaJonesCascade* pC;
  SViolaJonesPrimitive* pP;
  MIA_RESULT_CODE res=ERR_OK;
  float Weight;

    // check arguments
    if ((ppS==NULL)||(nama==NULL))
      return ERR_GEN_NULLPOINTER;
    if (*ppS!=NULL)
      return ERR_GEN_INITIALISED_ALREADY;
    // load file to buffer
    if ((res = FILE_AllocLoadFromFileExt((uint8**)&buf,(uint32*)&len,nama))!=ERR_OK)
      goto VJ_Load_exit;
    // add zero to the end
    if ((buf = (char*)realloc(buf,len+1))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto VJ_Load_exit;
    }
    buf[len] = 0;
    // allocate VJ
    if ((pS = (SViolaJonesClassifier*)malloc(sizeof(*pS)))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto VJ_Load_exit;
    }
    memset(pS,0,sizeof(*pS));
    // read common sizes
    if ((res = sscanf(buf+idx,"%d %d %d\n",&(pS->xs),&(pS->ys),&(pS->nCascade)))!=3)
    {
      res = ERR_GEN_BADDATATYPE;
      goto VJ_Load_exit;
    }
    if ((pS->xs<=0)||(pS->xs>256)||(pS->ys<=0)||(pS->ys>256)||
        (pS->nCascade<=0)||(pS->nCascade>100))
    {
      res = ERR_GEN_BAD_DATA;
      goto VJ_Load_exit;
    }
    xs = pS->xs;
    ys = pS->ys;
    // skip line
    for (;(buf[idx]!=13)&&(buf[idx]!=10)&&(idx<len);idx++);
    for (;((buf[idx]==13)||(buf[idx]==10))&&(idx<len);idx++);
    // allocate cascades
    if ((pS->pCascades = (SViolaJonesCascade*)malloc(sizeof(*pC)*pS->nCascade))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto VJ_Load_exit;
    }
    pC = pS->pCascades;
    memset(pC,0,sizeof(*pC)*pS->nCascade);
    // read cascades
    for (cCascade=0;cCascade<pS->nCascade;cCascade++)
    {
      pC = pS->pCascades+cCascade;
      // number of primitives
      if ((res = sscanf(buf+idx,"%d\n",&(pC->nPrimitives)))!=1)
      {
        res = ERR_GEN_BADDATATYPE;
        goto VJ_Load_exit;
      }
      if ((pC->nPrimitives<=0)||(pC->nPrimitives>100))
      {
        res = ERR_GEN_BAD_DATA;
        goto VJ_Load_exit;
      }
      // skip line
      for (;(buf[idx]!=13)&&(buf[idx]!=10)&&(idx<len);idx++);
      for (;((buf[idx]==13)||(buf[idx]==10))&&(idx<len);idx++);
      if ((pC->pPrimitives = (SViolaJonesPrimitive*)malloc(sizeof(*pP)*pC->nPrimitives))==NULL)
      {
        res = ERR_GEN_NOMEMORY;
        goto VJ_Load_exit;
      }
      pP = pC->pPrimitives;
      memset(pP,0,sizeof(*pP)*pC->nPrimitives);
      // read primitives
      for (cPrimitive=0;cPrimitive<pC->nPrimitives;cPrimitive++)
      {
        pP = pC->pPrimitives+cPrimitive;
        if ((res = sscanf(buf+idx,"%f   %d %d   %d %d %d %d  %d %d %d %d\n",
              &Weight,&(pP->nType),&(pP->nThreshold),
              pP->idx+0,pP->idx+1,pP->idx+2,pP->idx+3,
              pP->idx+4,pP->idx+5,pP->idx+6,pP->idx+7))!=11)
        {
          res = ERR_GEN_BADDATATYPE;
          goto VJ_Load_exit;
        }
        pP->weight = Weight;
        if (/*(pP->weight==0.)||*/(pP->nType<0)||(pP->nType>8)||
            (pP->idx[0]<0)||(pP->idx[0]>=(xs+1)*(ys+1))||
            (pP->idx[1]<0)||(pP->idx[1]>=(xs+1)*(ys+1))||
            (pP->idx[2]<0)||(pP->idx[2]>=(xs+1)*(ys+1))||
            (pP->idx[3]<0)||(pP->idx[3]>=(xs+1)*(ys+1))||
            (pP->idx[4]<0)||(pP->idx[4]>=(xs+1)*(ys+1))||
            (pP->idx[5]<0)||(pP->idx[5]>=(xs+1)*(ys+1))||
            (pP->idx[6]<0)||(pP->idx[6]>=(xs+1)*(ys+1))||
            (pP->idx[7]<0)||(pP->idx[7]>=(xs+1)*(ys+1)) )
        {
          res = ERR_GEN_BAD_DATA;
          goto VJ_Load_exit;
        }
        // skip line
        for (;(buf[idx]!=13)&&(buf[idx]!=10)&&(idx<len);idx++);
        for (;((buf[idx]==13)||(buf[idx]==10))&&(idx<len);idx++);
      }
    }
    res = ERR_OK;
VJ_Load_exit:;
    if (buf)
      free(buf);
    if (res!=ERR_OK)
      VJ_Free(&pS);
    *ppS = pS;
    return res;
}

MIA_RESULT_CODE VJ_Free(
  SViolaJonesClassifier** ppS)
{
  int cCascade;
  SViolaJonesClassifier* pS;
  SViolaJonesCascade* pC;

    if (ppS==NULL)
      return ERR_GEN_NULLPOINTER;
    pS = *ppS;
    if (pS==NULL)
      return ERR_OK;
    if (pS->pCascades)
    {
      for (cCascade=0;cCascade<pS->nCascade;cCascade++)
      {
        pC = pS->pCascades+cCascade;
        if (pC)
        {
//          if (pC->pWeights)
  //          free(pC->pWeights);
          if (pC->pPrimitives)
            free(pC->pPrimitives);
        }
      }
      free(pS->pCascades);
    }
    free(pS);
    *ppS = NULL;
    return ERR_OK;
}

MIA_RESULT_CODE VJ_Classify(
  int* pDecision,
  const SViolaJonesClassifier* pSVJ,
  const uint8* im,
  int str)
{
  int *iim=NULL,xs,ys,cPrimitive,cCascade,tmpint,deci=-1;
  SViolaJonesCascade *cas;
  SViolaJonesPrimitive *pri;
  double dValPos,dBestVal;
  MIA_RESULT_CODE res = ERR_OK;

    xs = pSVJ->xs;
    ys = pSVJ->ys;
    iim = (int*)malloc((xs+1)*(ys+1)*sizeof(*iim));
    IPL_PMP_MakeIntegralImage_prim(
      iim,  // OUT: integer image
      im,  // IN:  source image
      xs,ys);                   // IN:  size
    for (cCascade=0;cCascade<pSVJ->nCascade;cCascade++)
    {
      cas = pSVJ->pCascades+cCascade;
      dValPos = dBestVal = 0.;
      //++ this is a strong classifier
      for (cPrimitive=0;cPrimitive<cas->nPrimitives;cPrimitive++)
      {
        pri = cas->pPrimitives+cPrimitive;
        if ((pri->nType&(VJ_PRI_3SIDED|VJ_PRI_2SIDED))==VJ_PRI_3SIDED)
          tmpint = +iim[pri->idx[0]]-2*iim[pri->idx[1]]+2*iim[pri->idx[2]]-iim[pri->idx[3]]
                   -iim[pri->idx[4]]+2*iim[pri->idx[5]]-2*iim[pri->idx[6]]+iim[pri->idx[7]];
        else
          tmpint = +iim[pri->idx[0]]-2*iim[pri->idx[1]]+iim[pri->idx[2]]
                   -iim[pri->idx[3]]+2*iim[pri->idx[4]]-iim[pri->idx[5]];
tmpint -= THM*pri->nThreshold;
        if ( ((tmpint>0)&&(!(pri->nType&VJ_PRI_NEGATIVE))) ||
             ((tmpint<0)&&( (pri->nType&VJ_PRI_NEGATIVE))) ) 
          dBestVal += pri->weight;
        dValPos += pri->weight;
      }
      if (dBestVal*2<=dValPos)
      {
        deci = 0;
        goto VJ_Classify_exit;
      }
    }
    deci = 1;
VJ_Classify_exit:;
    if (iim)
      free(iim);
    *pDecision = deci;
    return res;
}
