/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  8 April 2013                                           */
/*      Revision: 1.0.00                                                 */
/*      Purpose:  linear chunks, merged by forest fire method            */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/
#define __FILENUM__ 118 // __FILENUM__TAG118

#include <stdio.h>
#include <malloc.h>
#include <string.h>
//#include "dbg.h"
#include "bpl.h"
#include "../ownipl.h"

#define NextColor(ko) (ko+2)
#define COCOBEG NextColor(0)
#define LINQTHR .2  // quality threshold of segment

#define AddPointToQueue_v1(dp) \
          if (cocoim[q+dp]==COCOBEG) \
          { \
            cocoim[queue[qinp] = q+dp] = ko; \
            qinp = (qinp+1)%QSZ; \
          }

#define AddPointToQueue_v2(dp) \
          if ((cocoim[q+dp]&1)&&(cocoim[q+dp]!=ko+1)) \
          { \
            cocoim[queue[qinp] = q+dp] = ko; \
            qinp = (qinp+1)%QSZ; \
          }

#define AddPointToQueue_v3(dp) \
          if ((cocoim[q+dp]&1)&&(cocoim[q+dp]!=ko+1)) \
          { \
            cocoim[queue[qinp] = q+dp] = ko; \
            qinp = (qinp+1)%QSZ; \
            nNoYieldRunlen = 0; \
          }

// get (a) linear (straight) chunks by (b) forest fire method by (c) coordinates
// if poor point is met, stop object enumeration by forest fire instantly
MIA_RESULT_CODE IPL_CHUNK_Linear_ForestFire_Coords_v1(
  CHUNK_Linear** ppsCL,     // OUT: array of segments
  int* pnCL,                // OUT: number of segments
  int32** ppnIdxs,          // OUT: list of point indexes (len = total number of enumerated points)
  int32** ppnCnts,          // OUT: list of index counts (len = number of enumerated objects+1)
  uint32** ppCoIm,          // OUT: 'colored' image of segments (32bpP)
  const unsigned char* im,  // IN:  source image
  const short* gxim,        // IN:  brightness gradients of source image
  const short* gyim,        //
  int xs,                   // IN:  source image size
  int ys,                   //
  int nMassThr)             // IN:  object mass threshold
{
  MIA_RESULT_CODE res = ERR_OK;
  int x,y,*queue=NULL,qinp,qout,p,q,xb,yb,nCL_u=0,nCL_a=0,QSZ,nPoints;
  int32 S1,Sx,Sy,Sxx,Sxy,Syy,Sgx,Sgy;
  double phi=0.,rho=0.,co,si,err1=0.,err2=0.;
  uint32 *cocoim=NULL,ko=COCOBEG;
  int32 *pnIdxs = NULL,*pnCnts = NULL;
  CHUNK_Linear* psCL = NULL,*pcl;

    // check arguments
    if ((ppsCL==NULL)||(pnCL==NULL)||(im==NULL)||(gxim==NULL)||(gyim==NULL)||
        (ppnIdxs==NULL)||(ppnCnts==NULL)||(ppCoIm==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((*ppsCL!=NULL)||(*ppCoIm!=NULL)||(*ppnCnts!=NULL)||(*ppnIdxs!=NULL))
      return ERR_GEN_INITIALISED_ALREADY;
    // init
    QSZ = xs+ys;
    // allocate
    if ((cocoim = (uint32*)malloc(sizeof(cocoim[0])*xs*ys))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto Linear_ForestFire_Coords_v1_exit;
    }
    if ((queue = (int*)malloc(QSZ*sizeof(queue[0])))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto Linear_ForestFire_Coords_v1_exit;
    }
    // fill color image: zero borders, initial value for objects
    nPoints = 0;
    memset(cocoim,0,xs*sizeof(cocoim[0]));
    memset(cocoim+(ys-1)*xs,0,xs*sizeof(cocoim[0]));
    for (y=1;y<ys-1;y++)
    {
      cocoim[y*xs] = cocoim[y*xs+xs-1] = 0;
      for (x=1;x<xs-1;x++)
        if (im[y*xs+x])
        {
          cocoim[y*xs+x] = COCOBEG;
          nPoints++;
        }
        else
          cocoim[y*xs+x] = 0;
    }
    // allocate arrays
    if ((pnIdxs = (int*)malloc(nPoints*sizeof(pnIdxs[0])))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto Linear_ForestFire_Coords_v1_exit;
    }
    if ((pnCnts = (int*)malloc((nPoints+1)*sizeof(pnCnts[0])))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto Linear_ForestFire_Coords_v1_exit;
    }
    // cycle in image
    nPoints = 0;
    for (p=xs*ys-1;p>=0;p--)
      if (cocoim[p]==COCOBEG)
      { // new object found
        ko = NextColor(ko);// kokoko
        qinp = qout = 0;
        q = p;
        pnCnts[nCL_u] = nPoints;
        // add initial point to the queue and delete from image
        AddPointToQueue_v1(0)
        // init the counters
        xb = q%xs;
        yb = q/xs;
        S1 = Sx = Sy = Sxx = Sxy = Syy = Sgx = Sgy = 0;
        // cycle by queue
        while ((qinp!=qout)&&((qout-qinp+QSZ)%QSZ>8))
        {
          // read point from queue
          q = queue[qout];
          qout = (qout+1)%QSZ;
          // get point coords
          x = q%xs-xb;
          y = q/xs-yb;
          // update counters
          S1++;
          Sx += x;
          Sy += y;
          Sxx += x*x;
          Sxy += x*y;
          Syy += y*y;
          Sgx += gxim[q];
          Sgy += gyim[q];
          if (S1>=3)
          { // calculate quality
            phi = .5*atan2(2*(S1*(double)Sxy-Sx*(double)Sy),
                           S1*((double)Sxx-Syy)-Sx*(double)Sx+Sy*(double)Sy);
            co = cos(phi);
            si = sin(phi);
            rho = (co*Sx+si*Sy)/S1;
            err1 = co*co*Sxx+2*co*si*Sxy+si*si*Syy-rho*rho*S1;
            rho = (-si*Sx+co*Sy)/S1;
            err2 = si*si*Sxx-2*co*si*Sxy+co*co*Syy-rho*rho*S1;
            if ((err1/S1>=LINQTHR)&&(err2/S1>=LINQTHR))
            { // last point makes object bad, undo the point
              // undo counters
              S1--;
              Sx -= x;
              Sy -= y;
              Sxx -= x*x;
              Sxy -= x*y;
              Syy -= y*y;
              Sgx -= gxim[q];
              Sgy -= gyim[q];
              // restore color of the point
              cocoim[q] = COCOBEG;
              // empty the queue, resetting color for points
              while (qinp!=qout)
              {
                q = queue[qout];
                qout = (qout+1)%QSZ;
                cocoim[q] = COCOBEG;
              }
              // terminate the queue cycle
              break;
            }
          }
          // point is good. add neighbors
          AddPointToQueue_v1(-xs)
          AddPointToQueue_v1( -1)
          AddPointToQueue_v1(  1)
          AddPointToQueue_v1( xs)
          AddPointToQueue_v1(-1-xs)
          AddPointToQueue_v1( 1-xs)
          AddPointToQueue_v1(-1+xs)
          AddPointToQueue_v1( 1+xs)
          // add point to list
          pnIdxs[nPoints++] = q;
        }
        if (qinp!=qout)
          break;  // queue overflow
        if (S1<nMassThr)
        { // small object - drop it
          nPoints = pnCnts[nCL_u];
          continue;
        }
        // save object data
        if (nCL_u==nCL_a)
          if ((psCL = (CHUNK_Linear*)realloc(psCL,(nCL_a+=1024)*sizeof(psCL[0])))==NULL)
          {
            res = ERR_GEN_NOMEMORY;
            goto Linear_ForestFire_Coords_v1_exit;
          }
        pcl = psCL+nCL_u;  // shortcut
        pcl->S1 = S1;
        pcl->id = ko;
        pcl->xc = xb+((double)Sx)/S1;
        pcl->yc = yb+((double)Sy)/S1;
        //pcl->gx = (Sgx+S1/2)/S1;
        //pcl->gy = (Sgy+S1/2)/S1;
        //pcl->c = ko;
        pcl->phi = phi;
        pcl->len = sqrt((double)S1);
        pcl->err = ((err1<err2)?err1:err2)/S1;
        nCL_u++;
      }
    if (p>=0)
    {
      res = ERR_GEN_INSUFFICIENT_BUFFER; // queue overflow
      goto Linear_ForestFire_Coords_v1_exit;
    }
    // set last (number of nCL_u+1) counter with total number of enumed points
    pnCnts[nCL_u] = nPoints;
Linear_ForestFire_Coords_v1_exit:;
    if (res!=ERR_OK)
    {
      if (psCL)
        free(psCL);
      psCL = NULL;
      nCL_u = 0;
      if (cocoim)
        free(cocoim);
      cocoim = NULL;
      if (pnIdxs)
        free(pnIdxs);
      pnIdxs = NULL;
      if (pnCnts)
        free(pnCnts);
      pnCnts = NULL;
    }
    else
    { // realloc to save space
      if ((psCL = (CHUNK_Linear*)realloc(psCL,sizeof(psCL[0])*nCL_u))==NULL)
      {
        res = ERR_GEN_NOMEMORY;
        goto Linear_ForestFire_Coords_v1_exit;
      }
      if ((pnCnts = (int32*)realloc(pnCnts,sizeof(pnCnts[0])*(nCL_u+1)))==NULL)
      {
        res = ERR_GEN_NOMEMORY;
        goto Linear_ForestFire_Coords_v1_exit;
      }
    }
    // release
    if (queue)
      free(queue);
    // output
    *ppsCL = psCL;
    *pnCL = nCL_u;
    *ppCoIm = cocoim;
    *ppnIdxs = pnIdxs;
    *ppnCnts = pnCnts;
    return res;
}

// get (a) linear (straight) chunks by (b) forest fire method by (c) coordinates
// if poor point is met, only this point is omitted
MIA_RESULT_CODE IPL_CHUNK_Linear_ForestFire_Coords_v2(
  CHUNK_Linear** ppsCL,     // OUT: array of segments
  int* pnCL,                // OUT: number of segments
  int32** ppnIdxs,          // OUT: list of point indexes (len = total number of enumerated points)
  int32** ppnCnts,          // OUT: list of index counts (len = number of enumerated objects+1)
  uint32** ppCoIm,          // OUT: 'colored' image of segments (32bpP)
  const unsigned char* im,  // IN:  source image
  const short* gxim,        // IN:  brightness gradients of source image
  const short* gyim,        //
  int xs,                   // IN:  source image size
  int ys,                   //
  int nMassThr)             // IN:  object mass threshold
{
  MIA_RESULT_CODE res = ERR_OK;
  int x,y,*queue=NULL,qinp,qout,p,q,xb,yb,nCL_u=0,nCL_a=0,QSZ,nPoints;
  int32 S1,Sx,Sy,Sxx,Sxy,Syy,Sgx,Sgy;
  double phi=0.,rho=0.,co,si,err1=0.,err2=0.;
  uint32 *cocoim=NULL,ko=COCOBEG;
  int32 *pnIdxs = NULL,*pnCnts = NULL;
  CHUNK_Linear* psCL = NULL,*pcl;

    // check arguments
    if ((ppsCL==NULL)||(pnCL==NULL)||(im==NULL)||(gxim==NULL)||(gyim==NULL)||
        (ppnIdxs==NULL)||(ppnCnts==NULL)||(ppCoIm==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((*ppsCL!=NULL)||(*ppCoIm!=NULL)||(*ppnCnts!=NULL)||(*ppnIdxs!=NULL))
      return ERR_GEN_INITIALISED_ALREADY;
    // init
    QSZ = xs+ys;
    // allocate
    if ((cocoim = (uint32*)malloc(sizeof(cocoim[0])*xs*ys))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto Linear_ForestFire_Coords_v2_exit;
    }
    if ((queue = (int*)malloc(QSZ*sizeof(queue[0])))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto Linear_ForestFire_Coords_v2_exit;
    }
    // fill color image: zero borders, initial value for objects
    nPoints = 0;
    memset(cocoim,0,xs*sizeof(cocoim[0]));
    memset(cocoim+(ys-1)*xs,0,xs*sizeof(cocoim[0]));
    for (y=1;y<ys-1;y++)
    {
      cocoim[y*xs] = cocoim[y*xs+xs-1] = 0;
      for (x=1;x<xs-1;x++)
        if (im[y*xs+x])
        {
          cocoim[y*xs+x] = 1;//COCOBEG;
          nPoints++;
        }
        else
          cocoim[y*xs+x] = 0;
    }
    // allocate arrays
    if ((pnIdxs = (int*)malloc(nPoints*sizeof(pnIdxs[0])))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto Linear_ForestFire_Coords_v2_exit;
    }
    if ((pnCnts = (int*)malloc((nPoints+1)*sizeof(pnCnts[0])))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto Linear_ForestFire_Coords_v2_exit;
    }
    // cycle in image
    nPoints = 0;
    for (p=xs*ys-1;p>=0;p--)
      if (cocoim[p]&1)//if (cocoim[p]==COCOBEG)
      { // new object found
        ko = NextColor(ko);// kokoko
        qinp = qout = 0;
        q = p;
        pnCnts[nCL_u] = nPoints;
        // add initial point to the queue and delete from image
        AddPointToQueue_v2(0)
        // init the counters
        xb = q%xs;
        yb = q/xs;
        S1 = Sx = Sy = Sxx = Sxy = Syy = Sgx = Sgy = 0;
        // cycle by queue
        while ((qinp!=qout)&&((qout-qinp+QSZ)%QSZ>8))
        {
          // read point from queue
          q = queue[qout];
          qout = (qout+1)%QSZ;
          // get point coords
          x = q%xs-xb;
          y = q/xs-yb;
          // update counters
          S1++;
          Sx += x;
          Sy += y;
          Sxx += x*x;
          Sxy += x*y;
          Syy += y*y;
          Sgx += gxim[q];
          Sgy += gyim[q];
          if (S1>=3)
          { // calculate quality
            phi = .5*atan2(2*(S1*(double)Sxy-Sx*(double)Sy),
                           S1*((double)Sxx-Syy)-Sx*(double)Sx+Sy*(double)Sy);
            co = cos(phi);
            si = sin(phi);
            rho = (co*Sx+si*Sy)/S1;
            err1 = co*co*Sxx+2*co*si*Sxy+si*si*Syy-rho*rho*S1;
            rho = (-si*Sx+co*Sy)/S1;
            err2 = si*si*Sxx-2*co*si*Sxy+co*co*Syy-rho*rho*S1;
            if ((err1/S1>=LINQTHR)&&(err2/S1>=LINQTHR))
            { // last point makes object bad, undo the point
              // undo counters
              S1--;
              Sx -= x;
              Sy -= y;
              Sxx -= x*x;
              Sxy -= x*y;
              Syy -= y*y;
              Sgx -= gxim[q];
              Sgy -= gyim[q];
              // restore color of the point
              cocoim[q]++;//cocoim[q] = COCOBEG;
              // proceed to next point
              continue;
            }
          }
          // point is good. add neighbors
          AddPointToQueue_v2(-xs)
          AddPointToQueue_v2( -1)
          AddPointToQueue_v2(  1)
          AddPointToQueue_v2( xs)
          AddPointToQueue_v2(-1-xs)
          AddPointToQueue_v2( 1-xs)
          AddPointToQueue_v2(-1+xs)
          AddPointToQueue_v2( 1+xs)
          // add point to list
          pnIdxs[nPoints++] = q;
        }
        if (qinp!=qout)
          break;  // queue overflow
        if (S1<nMassThr)
        { // small object - drop it
          nPoints = pnCnts[nCL_u];
          continue;
        }
        // save object data
        if (nCL_u==nCL_a)
          if ((psCL = (CHUNK_Linear*)realloc(psCL,(nCL_a+=1024)*sizeof(psCL[0])))==NULL)
          {
            res = ERR_GEN_NOMEMORY;
            goto Linear_ForestFire_Coords_v2_exit;
          }
        pcl = psCL+nCL_u;  // shortcut
        pcl->S1 = S1;
        pcl->id = ko;
        pcl->xc = xb+((double)Sx)/S1;
        pcl->yc = yb+((double)Sy)/S1;
        //pcl->gx = (Sgx+S1/2)/S1;
        //pcl->gy = (Sgy+S1/2)/S1;
        //pcl->c = ko;
        pcl->phi = phi;
        pcl->len = sqrt((double)S1);
        pcl->err = ((err1<err2)?err1:err2)/S1;
        nCL_u++;
      }
    if (p>=0)
    {
      res = ERR_GEN_INSUFFICIENT_BUFFER; // queue overflow
      goto Linear_ForestFire_Coords_v2_exit;
    }
    // set last (number of nCL_u+1) counter with total number of enumed points
    pnCnts[nCL_u] = nPoints;
Linear_ForestFire_Coords_v2_exit:;
    if (res!=ERR_OK)
    {
      if (psCL)
        free(psCL);
      psCL = NULL;
      nCL_u = 0;
      if (cocoim)
        free(cocoim);
      cocoim = NULL;
      if (pnIdxs)
        free(pnIdxs);
      pnIdxs = NULL;
      if (pnCnts)
        free(pnCnts);
      pnCnts = NULL;
    }
    else
    { // realloc to save space
      if ((psCL = (CHUNK_Linear*)realloc(psCL,sizeof(psCL[0])*nCL_u))==NULL)
      {
        res = ERR_GEN_NOMEMORY;
        goto Linear_ForestFire_Coords_v2_exit;
      }
      if ((pnCnts = (int32*)realloc(pnCnts,sizeof(pnCnts[0])*(nCL_u+1)))==NULL)
      {
        res = ERR_GEN_NOMEMORY;
        goto Linear_ForestFire_Coords_v2_exit;
      }
    }
    // release
    if (queue)
      free(queue);
    // output
    *ppsCL = psCL;
    *pnCL = nCL_u;
    *ppCoIm = cocoim;
    *ppnIdxs = pnIdxs;
    *ppnCnts = pnCnts;
    return res;
}

// get (a) linear (straight) chunks by (b) forest fire method by (c) coordinates
// if poor point is met, it is postponed to end of queue
MIA_RESULT_CODE IPL_CHUNK_Linear_ForestFire_Coords_v3(
  CHUNK_Linear** ppsCL,     // OUT: array of segments
  int* pnCL,                // OUT: number of segments
  int32** ppnIdxs,          // OUT: list of point indexes (len = total number of enumerated points)
  int32** ppnCnts,          // OUT: list of index counts (len = number of enumerated objects+1)
  uint32** ppCoIm,          // OUT: 'colored' image of segments (32bpP)
  const unsigned char* im,  // IN:  source image
  const short* gxim,        // IN:  brightness gradients of source image
  const short* gyim,        //
  int xs,                   // IN:  source image size
  int ys,                   //
  int nMassThr)             // IN:  object mass threshold
{
  MIA_RESULT_CODE res = ERR_OK;
  int x,y,*queue=NULL,qinp,qout,p,q,xb,yb,nCL_u=0,nCL_a=0,QSZ,nPoints;
  int nNoYieldRunlen=0;
  int32 S1,Sx,Sy,Sxx,Sxy,Syy,Sgx,Sgy;
  double phi=0.,rho=0.,co,si,err1=0.,err2=0.;
  uint32 *cocoim=NULL,ko=COCOBEG;
  int32 *pnIdxs = NULL,*pnCnts = NULL;
  CHUNK_Linear* psCL = NULL,*pcl;

    // check arguments
    if ((ppsCL==NULL)||(pnCL==NULL)||(im==NULL)||(gxim==NULL)||(gyim==NULL)||
        (ppnIdxs==NULL)||(ppnCnts==NULL)||(ppCoIm==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((*ppsCL!=NULL)||(*ppCoIm!=NULL)||(*ppnCnts!=NULL)||(*ppnIdxs!=NULL))
      return ERR_GEN_INITIALISED_ALREADY;
    // init
    QSZ = xs+ys;
    // allocate
    if ((cocoim = (uint32*)malloc(sizeof(cocoim[0])*xs*ys))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto Linear_ForestFire_Coords_v3_exit;
    }
    if ((queue = (int*)malloc(QSZ*sizeof(queue[0])))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto Linear_ForestFire_Coords_v3_exit;
    }
    // fill color image: zero borders, initial value for objects
    nPoints = 0;
    memset(cocoim,0,xs*sizeof(cocoim[0]));
    memset(cocoim+(ys-1)*xs,0,xs*sizeof(cocoim[0]));
    for (y=1;y<ys-1;y++)
    {
      cocoim[y*xs] = cocoim[y*xs+xs-1] = 0;
      for (x=1;x<xs-1;x++)
        if (im[y*xs+x])
        {
          cocoim[y*xs+x] = 1;//COCOBEG;
          nPoints++;
        }
        else
          cocoim[y*xs+x] = 0;
    }
    // allocate arrays
    if ((pnIdxs = (int*)malloc(nPoints*sizeof(pnIdxs[0])))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto Linear_ForestFire_Coords_v3_exit;
    }
    if ((pnCnts = (int*)malloc((nPoints+1)*sizeof(pnCnts[0])))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto Linear_ForestFire_Coords_v3_exit;
    }
    // cycle in image
    nPoints = 0;
    for (p=xs*ys-1;p>=0;p--)
      if (cocoim[p]&1)//if (cocoim[p]==COCOBEG)
      { // new object found
        ko = NextColor(ko);// kokoko
        qinp = qout = 0;
        q = p;
        pnCnts[nCL_u] = nPoints;
        // add initial point to the queue and delete from image
        AddPointToQueue_v3(0)
        // init the counters
        xb = q%xs;
        yb = q/xs;
        S1 = Sx = Sy = Sxx = Sxy = Syy = Sgx = Sgy = 0;
        // cycle by queue
        while ((qinp!=qout)&&((qout-qinp+QSZ)%QSZ>8)&&(nNoYieldRunlen<(qinp-qout+QSZ)%QSZ))
        {
          // read point from queue
          q = queue[qout];
          qout = (qout+1)%QSZ;
          // get point coords
          x = q%xs-xb;
          y = q/xs-yb;
          // update counters
          S1++;
          Sx += x;
          Sy += y;
          Sxx += x*x;
          Sxy += x*y;
          Syy += y*y;
          Sgx += gxim[q];
          Sgy += gyim[q];
          if (S1>=3)
          { // calculate quality
            phi = .5*atan2(2*(S1*(double)Sxy-Sx*(double)Sy),
                           S1*((double)Sxx-Syy)-Sx*(double)Sx+Sy*(double)Sy);
            co = cos(phi);
            si = sin(phi);
            rho = (co*Sx+si*Sy)/S1;
            err1 = co*co*Sxx+2*co*si*Sxy+si*si*Syy-rho*rho*S1;
            rho = (-si*Sx+co*Sy)/S1;
            err2 = si*si*Sxx-2*co*si*Sxy+co*co*Syy-rho*rho*S1;
            if ((err1/S1>=LINQTHR)&&(err2/S1>=LINQTHR))
            { // last point makes object bad, undo the point
              // undo counters
              S1--;
              Sx -= x;
              Sy -= y;
              Sxx -= x*x;
              Sxy -= x*y;
              Syy -= y*y;
              Sgx -= gxim[q];
              Sgy -= gyim[q];
              // postpone the point to end of queue
              queue[qinp] = q;
              qinp = (qinp+1)%QSZ;
              nNoYieldRunlen++;
              // proceed to next point from queue
              continue;
            }
          }
          // point is good. add neighbors
          AddPointToQueue_v3(-xs)
          AddPointToQueue_v3( -1)
          AddPointToQueue_v3(  1)
          AddPointToQueue_v3( xs)
          AddPointToQueue_v3(-1-xs)
          AddPointToQueue_v3( 1-xs)
          AddPointToQueue_v3(-1+xs)
          AddPointToQueue_v3( 1+xs)
          // add point to list
          pnIdxs[nPoints++] = q;
        }
        if ((qinp!=qout)&& (nNoYieldRunlen<(qinp-qout+QSZ)%QSZ)) 
          break;  // queue overflow
        if (!(nNoYieldRunlen<(qinp-qout+QSZ)%QSZ))
        { // some points did not merge, restore their initial color
          while (qinp!=qout)
          {
            q = queue[qout];
            qout = (qout+1)%QSZ;
            cocoim[q] = 1;
          }
        }
        if (S1<nMassThr)
        { // small object - drop it
          nPoints = pnCnts[nCL_u];
          continue;
        }
        // save object data
        if (nCL_u==nCL_a)
          if ((psCL = (CHUNK_Linear*)realloc(psCL,(nCL_a+=1024)*sizeof(psCL[0])))==NULL)
          {
            res = ERR_GEN_NOMEMORY;
            goto Linear_ForestFire_Coords_v3_exit;
          }
        pcl = psCL+nCL_u;  // shortcut
        pcl->S1 = S1;
        pcl->id = ko;
        pcl->xc = xb+((double)Sx)/S1;
        pcl->yc = yb+((double)Sy)/S1;
        //pcl->gx = (Sgx+S1/2)/S1;
        //pcl->gy = (Sgy+S1/2)/S1;
        //pcl->c = ko;
        pcl->phi = phi;
        pcl->len = sqrt((double)S1);
        pcl->err = ((err1<err2)?err1:err2)/S1;
        nCL_u++;
      }
    if (p>=0)
    {
      res = ERR_GEN_INSUFFICIENT_BUFFER; // queue overflow
      goto Linear_ForestFire_Coords_v3_exit;
    }
    // set last (number of nCL_u+1) counter with total number of enumed points
    pnCnts[nCL_u] = nPoints;
Linear_ForestFire_Coords_v3_exit:;
    if (res!=ERR_OK)
    {
      if (psCL)
        free(psCL);
      psCL = NULL;
      nCL_u = 0;
      if (cocoim)
        free(cocoim);
      cocoim = NULL;
      if (pnIdxs)
        free(pnIdxs);
      pnIdxs = NULL;
      if (pnCnts)
        free(pnCnts);
      pnCnts = NULL;
    }
    else
    { // realloc to save space
      if ((psCL = (CHUNK_Linear*)realloc(psCL,sizeof(psCL[0])*nCL_u))==NULL)
      {
        res = ERR_GEN_NOMEMORY;
        goto Linear_ForestFire_Coords_v3_exit;
      }
      if ((pnCnts = (int32*)realloc(pnCnts,sizeof(pnCnts[0])*(nCL_u+1)))==NULL)
      {
        res = ERR_GEN_NOMEMORY;
        goto Linear_ForestFire_Coords_v3_exit;
      }
    }
    // release
    if (queue)
      free(queue);
    // output
    *ppsCL = psCL;
    *pnCL = nCL_u;
    *ppCoIm = cocoim;
    *ppnIdxs = pnIdxs;
    *ppnCnts = pnCnts;
    return res;
}
