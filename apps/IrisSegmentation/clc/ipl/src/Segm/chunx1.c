/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  28 November 2012                                       */
/*      Revision: 1.0.00                                                 */
/*      Purpose:  construct chunks by points appending (forest fire)     */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/
#define __FILENUM__ 115 // __FILENUM__TAG115

#include <stdio.h>
#include <malloc.h>
#include <string.h>
//#include "dbg.h"
#include "bpl.h"
#include "../ownipl.h"

#define QSZ (1024*16)
#define kokoko ko = ko*1664525+1013904223;

#define COCOBEG 0xffff

#define ADDNB2(dp) \
          if (cocoim[q+dp]==COCOBEG) \
          { \
            cocoim[queue[qinp] = q+dp] = ko; \
            qinp = (qinp+1)%QSZ; \
          }

MIA_RESULT_CODE IPL_CONOB_ChunkDetect_v1(
  SChunkData** ppsCDO,      // OUT: array of segments
  int* pnCDO,               // OUT: number of segments
  const unsigned char* im,  // IN:  source image
  const short* gxim,        // IN:  brightness gradients of source image
  const short* gyim,        //
  int xs,                   // IN:  source image size
  int ys,                   //
  const char* nambeg)       // IN:  for debug, set to NULL
{
  MIA_RESULT_CODE res = ERR_OK;
  int x,y,*queue=NULL,qinp,qout,p,q,xb,yb,nCDO_u=0,nCDO_a=0;
  int32 S1,Sx,Sy,Sxx,Sxy,Syy,Sgx,Sgy;
  double phi,rho,co,si,err1,err2;//,LG;
  uint32* cocoim=NULL,ko=COCOBEG;
  SChunkData* psCDO = NULL,*pcdo;
  unsigned char *graim=NULL;
//  char nama[FILENAME_MAX];

    // check arguments
    if ((ppsCDO==NULL)||(pnCDO==NULL)||(im==NULL)||(gxim==NULL)||(gyim==NULL))
      return ERR_GEN_NULLPOINTER;
    if (*ppsCDO!=NULL)
      return ERR_GEN_INITIALISED_ALREADY;
    // allocate
    if ((cocoim = (uint32*)malloc(sizeof(cocoim[0])*xs*ys))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto ChuDex;
    }
    if ((queue = (int*)malloc(QSZ*sizeof(queue[0])))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto ChuDex;
    }
    if (nambeg)
    {
      graim = (unsigned char*)malloc(xs*ys*3);
      memset(graim,0,xs*ys*3);
    }
    // fill image
    memset(cocoim,0,xs*sizeof(cocoim[0]));
    memset(cocoim+(ys-1)*xs,0,xs*sizeof(cocoim[0]));
    for (y=1;y<ys-1;y++)
    {
      cocoim[y*xs] = cocoim[y*xs+xs-1] = 0;
      for (x=1;x<xs-1;x++)
        cocoim[y*xs+x] = (im[y*xs+x])?COCOBEG:0;
    }
    // cycle in image
    for (p=xs*ys-1;p>=0;p--)
      if (cocoim[p]==COCOBEG)
      { // new object found
        kokoko
        qinp = qout = 0;
        q = p;
        // add initial point to the queue and delete from image
        ADDNB2(0)
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
            if ((err1/S1>=.2)&&(err2/S1>=.2))
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
              cocoim[q] = COCOBEG;
              // empty the queue, resetting point color
              while (qinp!=qout)
              {
                q = queue[qout];
                qout = (qout+1)%QSZ;
                cocoim[q] = COCOBEG;
              }
              //continue; // proceed to next point from queue
              break;
            }
          }
          // point is good. add neighbors
          ADDNB2(-xs)
          ADDNB2( -1)
          ADDNB2(  1)
          ADDNB2( xs)
          ADDNB2(-1-xs)
          ADDNB2( 1-xs)
          ADDNB2(-1+xs)
          ADDNB2( 1+xs)
        }
        if (qinp!=qout)
          break;  // queue overflow
        if (S1<10)
          continue;
//ch printf("N,x,y,gx,gy=(%d,%d,%d,%d,%d)\n",
  //     S1,xb+(Sx+S1/2)/S1,yb+(Sy+S1/2)/S1,(Sgx+S1/2)/S1,(Sgy+S1/2)/S1);
        // save object data
        if (nCDO_u==nCDO_a)
          if ((psCDO = (SChunkData*)realloc(psCDO,(nCDO_a+=1024)*sizeof(psCDO[0])))==NULL)
          {
            res = ERR_GEN_NOMEMORY;
            goto ChuDex;
          }
        pcdo = psCDO+nCDO_u;  // shortcut
        pcdo->n = S1;
        pcdo->x = xb+(Sx+S1/2)/S1;
        pcdo->y = yb+(Sy+S1/2)/S1;
        pcdo->gx = (Sgx+S1/2)/S1;
        pcdo->gy = (Sgy+S1/2)/S1;
        pcdo->c = ko;
/*        if (nambeg)
        { // draw in image
          LG = sqrt(pcdo->gx*pcdo->gx+pcdo->gy*pcdo->gy);
          DBGL_DRAW_LineInRGB(graim,xs*3,xs,ys,pcdo->c,pcdo->x,pcdo->y,
                              pcdo->x+(int)(pcdo->gx*10/LG),
                              pcdo->y+(int)(pcdo->gy*10/LG));
          DBGL_DRAW_RectInRGB(graim,xs*3,xs,ys,pcdo->c,pcdo->x-1,pcdo->y-1,3,3);
        }*/
        nCDO_u++;
      }
    if (p>=0)
    {
      res = ERR_GEN_INSUFFICIENT_BUFFER; // queue overflow
      goto ChuDex;
    }
/*    if (nambeg)
    {
      // save image of segments
      sprintf(nama,"nCDO_u=%d",nCDO_u);
      DBGL_TYPE_TextInRGBA(cocoim,xs,ys,nama,0,0,0xff00,-1);
      strcpy(nama,nambeg);
      strcat(nama,"_1_segments.bmp");
      DBGL_FILE_SaveRGBA8Image(cocoim,xs,xs,ys,nama);
      // save image of gradients
      strcpy(nama,nambeg);
      strcat(nama,"_2_gradients.bmp");
      DBGL_FILE_SaveRGB8Image(graim,xs,xs,ys,nama);
    }*/
ChuDex:;
//ch printf("nCDO_u=%d\n",nCDO_u);
    if (res!=ERR_OK)
    {
      free(psCDO);
      psCDO = NULL;
      nCDO_u = 0;
    }
    // output
    *ppsCDO = psCDO;
    *pnCDO = nCDO_u;
    // release
    if (cocoim)
      free(cocoim);
    if (graim)
      free(graim);
    if (queue)
      free(queue);
    return res;
}

MIA_RESULT_CODE IPL_CONOB_ArchunkDetect_v1(
  SArchunkData** ppsCDI,    // OUT: array of arcs
  int* pnCDI,               // OUT: number of segments
  const unsigned char* im,  // IN:  source image
  const short* gxim,        // IN:  brightness gradients of source image
  const short* gyim,        //
  int xs,                   // IN:  source image size
  int ys,                   //
  const char* nambeg)       // IN:  for debug, set to NULL
{
  MIA_RESULT_CODE res = ERR_OK;
  int x,y,*queue=NULL,qinp,qout,p,q,xb,yb,nCDI_u=0,nCDI_a=0,isbreak;
  int64 S1,Sx,Sy,Sxx,Sxy,Syy,Sxxx,Sxxy,Sxyy,Syyy,Sxxxx,Sxxyy,Syyyy;
  double phi,rho,co,si,err1,err2;
  uint32* cocoim=NULL,ko=COCOBEG;
  SArchunkData* psCDI = NULL,*pcdi;
  unsigned char *graim=NULL;
//  char nama[FILENAME_MAX];
  SLinearSystem pLinSys3;
  double u,x0,y0;

    // check arguments
    if ((ppsCDI==NULL)||(pnCDI==NULL)||(im==NULL)||(gxim==NULL)||(gyim==NULL))
      return ERR_GEN_NULLPOINTER;
    if (*ppsCDI!=NULL)
      return ERR_GEN_INITIALISED_ALREADY;
    // allocate
    if ((cocoim = (uint32*)malloc(sizeof(cocoim[0])*xs*ys))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto ChuDex;
    }
    if ((queue = (int*)malloc(QSZ*sizeof(queue[0])))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto ChuDex;
    }
    if ((res = BPL_LIN_AllocateSystem(&pLinSys3,3,3))!=ERR_OK)
    {
      res = ERR_GEN_NOMEMORY;
      goto ChuDex;
    }
    if (nambeg)
    {
      graim = (unsigned char*)malloc(xs*ys*3);
      memset(graim,0,xs*ys*3);
/*      strcpy(nama,nambeg);
      strcat(nama,"_0_source.bmp");
      DBGL_FILE_SaveUint8Image(im,xs,xs,ys,nama);*/
    }
    // fill image
    memset(cocoim,0,xs*sizeof(cocoim[0]));
    memset(cocoim+(ys-1)*xs,0,xs*sizeof(cocoim[0]));
    for (y=1;y<ys-1;y++)
    {
      cocoim[y*xs] = cocoim[y*xs+xs-1] = 0;
      for (x=1;x<xs-1;x++)
        cocoim[y*xs+x] = (im[y*xs+x])?COCOBEG:0;
    }
    // cycle in image
    for (p=xs*ys-1;p>=0;p--)
      if (cocoim[p]==COCOBEG)
      { // new object found
        kokoko
        qinp = qout = 0;
        q = p;
        // add initial point to the queue and delete from image
        ADDNB2(0)
        // init the counters
        xb = q%xs;
        yb = q/xs;
        S1 = Sx = Sy = Sxx = Sxy = Syy = Sxxx = Sxxy = Sxyy = Syyy = Sxxxx = Sxxyy = Syyyy = 0;
        u = x0 = y0 = 0;
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
          Sxxx += x*x*x;
          Sxxy += x*x*y;
          Sxyy += x*y*y;
          Syyy += y*y*y;
          Sxxxx += x*x*x*x;
          Sxxyy += x*x*y*y;
          Syyyy += y*y*y*y;
          isbreak = 0;
          if (S1>=20)
          { // enough point to verify circularity
            pLinSys3.A[0][0] = (double)(2*S1);
            pLinSys3.A[0][1] = (double)(2*Sx);
            pLinSys3.A[0][2] = (double)(2*Sy);
            pLinSys3.A[1][0] = (double)(2*Sx);
            pLinSys3.A[1][1] = (double)(2*Sxx);
            pLinSys3.A[1][2] = (double)(2*Sxy);
            pLinSys3.A[2][0] = (double)(2*Sy);
            pLinSys3.A[2][1] = (double)(2*Sxy);
            pLinSys3.A[2][2] = (double)(2*Syy);
            pLinSys3.b[0] = (double)(Sxx+Syy);
            pLinSys3.b[1] = (double)(Sxxx+Sxyy);
            pLinSys3.b[2] = (double)(Sxxy+Syyy);
            BPL_LIN_SolveSystem(&pLinSys3);
            u  = pLinSys3.x[0];
            x0 = pLinSys3.x[1];
            y0 = pLinSys3.x[2];
            /*err1 =  S1*u*u+Sxx*x0*x0+Syy*y0*y0
                   +2*Sx*u*x0+2*Sy*u*y0+2*Sxy*x0*y0
                   -(Sxx+Syy)*u-(Sxxx+Sxyy)*x0-(Sxxy+Syyy)*y0
                   +.25*(Sxxxx+2*Sxxyy+Syyyy);
            if (err1/S1>10.)
              goto remove_point;*/
            err1 = (x-x0)*(x-x0)+(y-y0)*(y-y0);
            err1 = sqrt(2*u+x0*x0+y0*y0)-sqrt(err1);
            if ((err1>5.)||(err1<-5.))
              goto remove_point;
          }
          else
            if (S1>=3)
            { // enough points to verify straightness
              phi = .5*atan2(2*(S1*(double)Sxy-Sx*(double)Sy),
                             S1*((double)Sxx-Syy)-Sx*(double)Sx+Sy*(double)Sy);
              co = cos(phi);
              si = sin(phi);
              rho = (co*Sx+si*Sy)/S1;
              err1 = co*co*Sxx+2*co*si*Sxy+si*si*Syy-rho*rho*S1;
              rho = (-si*Sx+co*Sy)/S1;
              err2 = si*si*Sxx-2*co*si*Sxy+co*co*Syy-rho*rho*S1;
              if ((err1/S1>=1)&&(err2/S1>=1))
                // last point makes object bad, undo the point
                goto remove_point;
            }
          // point is good. add neighbors
          ADDNB2(-xs)
          ADDNB2( -1)
          ADDNB2(  1)
          ADDNB2( xs)
          ADDNB2(-1-xs)
          ADDNB2( 1-xs)
          ADDNB2(-1+xs)
          ADDNB2( 1+xs)
          continue;
remove_point:;
          // undo counters
          S1--;
          Sx -= x;
          Sy -= y;
          Sxx -= x*x;
          Sxy -= x*y;
          Syy -= y*y;
          Sxxx -= x*x*x;
          Sxxy -= x*x*y;
          Sxyy -= x*y*y;
          Syyy -= y*y*y;
          Sxxxx -= x*x*x*x;
          Sxxyy -= x*x*y*y;
          Syyyy -= y*y*y*y;
          cocoim[q] = COCOBEG;
          // empty the queue, resetting point color
          while (qinp!=qout)
          {
            q = queue[qout];
            qout = (qout+1)%QSZ;
            cocoim[q] = COCOBEG;
          }
          //continue; // proceed to next point from queue
          break;
        }
        if (qinp!=qout)
          break;  // queue overflow
        if (S1<20)
          continue;
        // save object data
        if (nCDI_u==nCDI_a)
          if ((psCDI = (SArchunkData*)realloc(psCDI,(nCDI_a+=1024)*sizeof(psCDI[0])))==NULL)
          {
            res = ERR_GEN_NOMEMORY;
            goto ChuDex;
          }
        pcdi = psCDI+nCDI_u;  // shortcut
        pcdi->n = (int)S1;
        pcdi->x = xb+(int)(x0+.5);
        pcdi->y = yb+(int)(y0+.5);
        pcdi->r = (int)(sqrt(2*u+x0*x0+y0*y0)+.5);
        pcdi->c = ko;
/*        if (nambeg)
        { // draw in image
          DBGL_DRAW_PlusInRGB24(graim,xs,ys,pcdi->x,pcdi->y,pcdi->r,pcdi->c);
        }*/
        nCDI_u++;
      }
    if (p>=0)
    {
      res = ERR_GEN_INSUFFICIENT_BUFFER; // queue overflow
      goto ChuDex;
    }
/*    if (nambeg)
    {
      // save image of segments
      sprintf(nama,"nCDI_u=%d",nCDI_u);
      DBGL_TYPE_TextInRGBA(cocoim,xs,ys,nama,0,0,0xff00,-1);
      strcpy(nama,nambeg);
      strcat(nama,"_1_segments.bmp");
      DBGL_FILE_SaveRGBA8Image(cocoim,xs,xs,ys,nama);
      // save image of gradients
//      strcpy(nama,nambeg);
  //    strcat(nama,"_2_gradients.bmp");
    //  DBGL_FILE_SaveRGB8Image(graim,xs,xs,ys,nama);
    }*/
ChuDex:;
//ch printf("nCDO_u=%d\n",nCDO_u);
    if (res!=ERR_OK)
    {
      free(psCDI);
      psCDI = NULL;
      nCDI_u = 0;
    }
    BPL_LIN_FreeSystem(&pLinSys3);
    // output
    *ppsCDI = psCDI;
    *pnCDI = nCDI_u;
    // release
    if (cocoim)
      free(cocoim);
    if (graim)
      free(graim);
    if (queue)
      free(queue);
    return res;
}
