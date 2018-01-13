/*------------------------------------------------------------------------*/
/*                                                                        */
/*      Created:  22 November 2013                                        */
/*      Revision: 1.0.00                                                  */
/*      Purpose:  EdgeDraw filter (inspired by Cihan Topal &              */
/*                                             Cuneyt Akinlar work)       */
/*      Authors:                                                          */
/*        Ivan Matveev                                                    */
/*      Also refer to http://ceng.anadolu.edu.tr/cv/EdgeDrawing/          */
/*------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dbg.h"
#include "bpl.h"
#include "stddefs.h"
#include "errorcodes.h"
#include "ipl.h"

#define V2 1
#define V1 4
#define V0 6
#define MULTI_HOR (V0+2*V1+2*V2)
#define MULTI_VER ((MULTI_HOR*(V0+2*V1+2*V2))/4)
#define GRATHR (MULTI_VER*3)
#define DIVG 4
#define VS ((V0)+2*(V1)+2*(V2))
#define SCALEGXY(p,n) ((p)-(n))
//#define LENGTHGXY(gx,gy) (short)(sqrt(((double)gx)*gx+((double)gy)*gy)+.5)
//#define LENGTHGXY(gx,gy) MIA_INT_sqrt((int)(gx)*(gx)+(int)(gy)*(gy))
//#define LENGTHGXY(gx,gy) (gx+gy)
#define LENGTHGXY(gx,gy) (gx+gy+ ((gx>2*gy)?(gx/3):((gy>2*gx)?(gy/3):0)) )
#define DIRGRAD(gx,gy,V,T) (((V)<(T))?0: ((gx)>(gy))?127:(-127) )
#define ANCHORTHR 2

//const int nMULTI_VER = MULTI_VER;

void IPL_PMP_SmoothX_short(
  short* sm1,
  const unsigned char* im,
  int xs,
  int ys)
{
  int y,idx,idx_fin;

    idx_fin = -2;
    for (y=0;y<ys;y++)
    {
      idx = idx_fin+4;
      idx_fin += xs;
      sm1[idx-2] = (V0+V1)*im[idx-2]+(V1+V2)*im[idx-1]+V2*im[idx];
      sm1[idx-1] = (V1+V2)*im[idx-2]+V0*im[idx-1]+V1*im[idx]+V2*im[idx+1];
      for (;idx<idx_fin;idx++)
        sm1[idx] = V2*im[idx-2]+V1*im[idx-1]+V0*im[idx]+V1*im[idx+1]+V2*im[idx+2];
      sm1[idx] = V2*im[idx-2]+V1*im[idx-1]+V0*im[idx]+(V1+V2)*im[idx+1];
      sm1[idx+1] = V2*im[idx-1]+(V1+V2)*im[idx]+(V0+V1)*im[idx+1];
    }
}

void IPL_PMP_SmoothY_short(
  short* sm1,
  const short* im,
  int xs,
  int ys)
{
  int x,idx,idx_fin;

    idx_fin = (ys-2)*xs;
    for (x=0;x<xs;x++)
    {
      idx = idx_fin-xs*(ys-4);
      sm1[idx-2*xs] = ((V0+V1)*im[idx-2*xs]+(V1+V2)*im[idx-1*xs]+V2*im[idx])/DIVG;
      sm1[idx-1*xs] = ((V1+V2)*im[idx-2*xs]+V0*im[idx-1*xs]+V1*im[idx]+V2*im[idx+1*xs])/DIVG;
      for (;idx<idx_fin;idx+=xs)
        sm1[idx] = (V2*im[idx-2*xs]+V1*im[idx-1*xs]+V0*im[idx]+V1*im[idx+1*xs]+V2*im[idx+2*xs])/DIVG;
      sm1[idx] = (V2*im[idx-2*xs]+V1*im[idx-xs]+V0*im[idx]+(V1+V2)*im[idx+xs])/DIVG;
      sm1[idx+xs] = (V2*im[idx-xs]+(V1+V2)*im[idx]+(V0+V1)*im[idx+xs])/DIVG;
      idx_fin++;
    }
}

void IPL_PMP_SmoothX_int(
  int* sm1,
  const int* im,
  int xs,
  int ys)
{
  int y,idx,idx_fin;

    idx_fin = -2;
    for (y=0;y<ys;y++)
    {
      idx = idx_fin+4;
      idx_fin += xs;
      sm1[idx-2] = (V0+V1)*im[idx-2]+(V1+V2)*im[idx-1]+V2*im[idx];
      sm1[idx-1] = (V1+V2)*im[idx-2]+V0*im[idx-1]+V1*im[idx]+V2*im[idx+1];
      for (;idx<idx_fin;idx++)
        sm1[idx] = V2*im[idx-2]+V1*im[idx-1]+V0*im[idx]+V1*im[idx+1]+V2*im[idx+2];
      sm1[idx] = V2*im[idx-2]+V1*im[idx-1]+V0*im[idx]+(V1+V2)*im[idx+1];
      sm1[idx+1] = V2*im[idx-1]+(V1+V2)*im[idx]+(V0+V1)*im[idx+1];
    }
}

void IPL_PMP_SmoothX_intchar(
  int* sm1,
  const uint8* im,
  int xs,
  int ys)
{
  int y,idx,idx_fin;

    idx_fin = -2;
    for (y=0;y<ys;y++)
    {
      idx = idx_fin+4;
      idx_fin += xs;
      sm1[idx-2] = (V0+V1)*im[idx-2]+(V1+V2)*im[idx-1]+V2*im[idx];
      sm1[idx-1] = (V1+V2)*im[idx-2]+V0*im[idx-1]+V1*im[idx]+V2*im[idx+1];
      for (;idx<idx_fin;idx++)
        sm1[idx] = V2*im[idx-2]+V1*im[idx-1]+V0*im[idx]+V1*im[idx+1]+V2*im[idx+2];
      sm1[idx] = V2*im[idx-2]+V1*im[idx-1]+V0*im[idx]+(V1+V2)*im[idx+1];
      sm1[idx+1] = V2*im[idx-1]+(V1+V2)*im[idx]+(V0+V1)*im[idx+1];
    }
}

void IPL_PMP_SmoothY_int(
  int* sm1,
  const int* im,
  int xs,
  int ys)
{
  int x,idx,idx_fin;

    idx_fin = (ys-2)*xs;
    for (x=0;x<xs;x++)
    {
      idx = idx_fin-xs*(ys-4);
      sm1[idx-2*xs] = ((V0+V1)*im[idx-2*xs]+(V1+V2)*im[idx-1*xs]+V2*im[idx])/DIVG;
      sm1[idx-1*xs] = ((V1+V2)*im[idx-2*xs]+V0*im[idx-1*xs]+V1*im[idx]+V2*im[idx+1*xs])/DIVG;
      for (;idx<idx_fin;idx+=xs)
        sm1[idx] = (V2*im[idx-2*xs]+V1*im[idx-1*xs]+V0*im[idx]+V1*im[idx+1*xs]+V2*im[idx+2*xs])/DIVG;
      sm1[idx] = (V2*im[idx-2*xs]+V1*im[idx-xs]+V0*im[idx]+(V1+V2)*im[idx+xs])/DIVG;
      sm1[idx+xs] = (V2*im[idx-xs]+(V1+V2)*im[idx]+(V0+V1)*im[idx+xs])/DIVG;
      idx_fin++;
    }
}

staticFunc void ownIPL_gracal(
  short* gx,
  short* gy,
  short* G,
  char* D,
  const short* im,
  int xs,
  int ys,
  int grathr)
{
  int y,idx,idx_fin,vx,vy;

    idx = 0;
    //== top line
    //-- left border
    gx[idx] = (short)(vx = SCALEGXY(im[idx+1],im[idx]));
    gy[idx] = (short)(vy = SCALEGXY(im[idx+xs],im[idx]));
    if (vx<0)
      vx = -vx;
    if (vy<0)
      vy = -vy;
    G[idx] = (short)(LENGTHGXY(vx,vy));
    D[idx] = DIRGRAD(vx,vy,G[idx],grathr);
    idx++;
    //-- line body
    idx_fin = idx+xs-2;
    for (;idx<idx_fin;idx++)
    {
      gx[idx] = (short)(vx = SCALEGXY(im[idx+1],im[idx-1]));
      gy[idx] = (short)(vy = SCALEGXY(im[idx+xs],im[idx]));
      if (vx<0)
        vx = -vx;
      if (vy<0)
        vy = -vy;
      G[idx] = (short)(LENGTHGXY(vx,vy));
      D[idx] = DIRGRAD(vx,vy,G[idx],grathr);
    }
    //--right border
    gx[idx] = (short)(vx = SCALEGXY(im[idx],im[idx-1]));
    gy[idx] = (short)(vy = SCALEGXY(im[idx+xs],im[idx]));
    if (vx<0)
      vx = -vx;
    if (vy<0)
      vy = -vy;
    G[idx] = (short)(LENGTHGXY(vx,vy));
    D[idx] = DIRGRAD(vx,vy,G[idx],grathr);
    idx++;
    //== line body
    for (y=1;y<ys-1;y++)
    {
      //-- left border
      gx[idx] = (short)(vx = SCALEGXY(im[idx+1],im[idx]));
      gy[idx] = (short)(vy = SCALEGXY(im[idx+xs],im[idx-xs]));
      if (vx<0)
        vx = -vx;
      if (vy<0)
        vy = -vy;
      G[idx] = (short)(LENGTHGXY(vx,vy));
      D[idx] = DIRGRAD(vx,vy,G[idx],grathr);
      idx++;
      //-- line body
      idx_fin = idx+xs-2;
      for (;idx<idx_fin;idx++)
      {
        gx[idx] = (short)(vx = SCALEGXY(im[idx+1],im[idx-1]));
        gy[idx] = (short)(vy = SCALEGXY(im[idx+xs],im[idx-xs]));
        if (vx<0)
          vx = -vx;
        if (vy<0)
          vy = -vy;
        G[idx] = (short)(LENGTHGXY(vx,vy));
        D[idx] = DIRGRAD(vx,vy,G[idx],grathr);
      }
      // right border
      gx[idx] = (short)(vx = SCALEGXY(im[idx],im[idx-1]));
      gy[idx] = (short)(vy = SCALEGXY(im[idx+xs],im[idx-xs]));
      if (vx<0)
        vx = -vx;
      if (vy<0)
        vy = -vy;
      G[idx] = (short)(LENGTHGXY(vx,vy));
      D[idx] = DIRGRAD(vx,vy,G[idx],grathr);
      idx++;
    }
    //== bottom line
    //-- left border
    gx[idx] = (short)(vx = SCALEGXY(im[idx+1],im[idx]));
    gy[idx] = (short)(vy = SCALEGXY(im[idx],im[idx-xs]));
    if (vx<0)
      vx = -vx;
    if (vy<0)
      vy = -vy;
    G[idx] = (short)(LENGTHGXY(vx,vy));
    D[idx] = DIRGRAD(vx,vy,G[idx],grathr);
    idx++;
    //-- line body
    idx_fin = idx+xs-2;
    for (;idx<idx_fin;idx++)
    {
      gx[idx] = (short)(vx = SCALEGXY(im[idx+1],im[idx-1]));
      gy[idx] = (short)(vy = SCALEGXY(im[idx],im[idx-xs]));
      if (vx<0)
        vx = -vx;
      if (vy<0)
        vy = -vy;
      G[idx] = (short)(LENGTHGXY(vx,vy));
      D[idx] = DIRGRAD(vx,vy,G[idx],grathr);
    }
    //--right border
    gx[idx] = (short)(vx = SCALEGXY(im[idx],im[idx-1]));
    gy[idx] = (short)(vy = SCALEGXY(im[idx],im[idx-xs]));
    if (vx<0)
      vx = -vx;
    if (vy<0)
      vy = -vy;
    G[idx] = (short)(LENGTHGXY(vx,vy));
    D[idx] = DIRGRAD(vx,vy,G[idx],grathr);
    idx++;
}

typedef struct
{
  int idx;
  int val;
} SAnchor;

staticFunc void ownIPL_getanchors(
  int* anchorlist,
  int* anchornum,
  char* anchorim,
  const short* G,
  const char* D,
  int xs,
  int ys,
  int dx,
  int dy)
{
  int y,idx,idx_fin,an=0;

    memset(anchorim,0,ys*xs);
    //== body
    for (y=(dy/2>1)?(dy/2):1;y<ys-1;y+=dy)
    {
      //-- line body
      idx_fin = y*xs+xs-1;
      for (idx=y*xs+((dx/2>1)?(dx/2):1);idx<idx_fin;idx+=dx)
      {
        //anchorim[idx] = 0;
        if (D[idx])
        {
          if (D[idx]<0)
          { // vertical, compare with vertical neigbs
            if ((G[idx]>G[idx-xs]+ANCHORTHR)&&(G[idx]>G[idx+xs]+ANCHORTHR))
            {
              anchorlist[an*2+0] = idx;
              anchorlist[an*2+1] = G[idx];
              anchorim[idx] = -127;  // vertical anchor
              an++;
            }
          }
          else
          { // horizontal, compare with horizontal neigbs
            if ((G[idx]>G[idx-1]+ANCHORTHR)&&(G[idx]>G[idx+1]+ANCHORTHR))
            {
              anchorlist[an*2+0] = idx;
              anchorlist[an*2+1] = G[idx];
              anchorim[idx] = 127;  // horizontal anchor
              an++;
            }
          }
        }
      }
    }
    // return
    *anchornum = an;
}
/*
staticFunc int sortanch(
  const void* e1,
  const void* e2)
{
    return ((int*)e2)[1]-((int*)e1)[1];
}
*/
#define cmp_anch_decr(anc1,anc2) ((anc1.val)>(anc2.val))
IMPLEMENT_HEAPSORT(ownIPL_heapSortAnchsDecreasing,SAnchor,cmp_anch_decr)
/*
staticFunc void trackrt(
  char* edge,
  int str,
  int idx,
  const short* G,
  const char* D,
  int isPos)
{
//    edge[idx] = 0;  // unmark edge
    edge[idx] = 127;
    // first step
    if (D[idx]>0) 
    { // horizontal edge (vertical gradient)
      if (isPos>=0)
      { // step right
        if ((G[idx+1]>=G[idx+str+1])&&(G[idx+1]>=G[idx-str+1]))
          idx += 1;  // =>right
        else
          if ((G[idx+str+1]>=G[idx+1])&&(G[idx+str+1]>=G[idx-str+1]))
            idx += str+1; // =>right-up
          else
            idx += -str+1;  // =>right-down
      }
      else
      { // step left
        if ((G[idx-1]>=G[idx+str-1])&&(G[idx-1]>=G[idx-str-1]))
          idx += -1;  // =>left
        else
          if ((G[idx+str-1]>=G[idx-1])&&(G[idx+str-1]>=G[idx-str-1]))
            idx += str-1; // =>left-up
          else
            idx += -str-1;  // =>left-down
      }
    }
    else
    { // vertical edge (horizontal gradient)
      if (isPos>=0)
      { // step up
        if ((G[idx+str]>=G[idx+1+str])&&(G[idx+str]>=G[idx-1+str]))
          idx += str;  // =>up
        else
          if ((G[idx+1+str]>=G[idx+str])&&(G[idx+1+str]>=G[idx-1+str]))
            idx += 1+str; // =>up-right
          else
            idx += -1+str;  // =>up-left
      }
      else
      { // step down
        if ((G[idx-str]>=G[idx+1-str])&&(G[idx-str]>=G[idx-1-str]))
          idx += -str;  // =>down
        else
          if ((G[idx+1-str]>=G[idx-str])&&(G[idx+1-str]>=G[idx-1-str]))
            idx += 1-str; // =>down-right
          else
            idx += -1-str;  // =>down-left
      }
    }
    // normal run now
    while ((!edge[idx])&&(D[idx]))
    {
      if (D[idx]>0) 
      { // horizontal edge (vertical gradient), walk horizontally
        { // step right
          if ((G[idx+1]>=G[idx+str+1])&&(G[idx+1]>=G[idx-str+1]))
            idx += 1;  // =>right
          else
            if ((G[idx+str+1]>=G[idx+1])&&(G[idx+str+1]>=G[idx-str+1]))
              idx += str+1; // =>right-up
            else
              idx += -str+1;  // =>right-down
        }
        else
        { // step left
          if ((G[idx-1]>=G[idx+str-1])&&(G[idx-1]>=G[idx-str-1]))
            idx += -1;  // =>left
          else
            if ((G[idx+str-1]>=G[idx-1])&&(G[idx+str-1]>=G[idx-str-1]))
              idx += str-1; // =>left-up
            else
              idx += -str-1;  // =>left-down
        }
      }
      else
      { // vertical edge (horizontal gradient), walk vertically
        if (isPos>=0)
        { // step up
          if ((G[idx+str]>=G[idx+1+str])&&(G[idx+str]>=G[idx-1+str]))
            idx += str;  // =>up
          else
            if ((G[idx+1+str]>=G[idx+str])&&(G[idx+1+str]>=G[idx-1+str]))
              idx += 1+str; // =>up-right
            else
              idx += -1+str;  // =>up-left
        }
        else
        { // step down
          if ((G[idx-str]>=G[idx+1-str])&&(G[idx-str]>=G[idx-1-str]))
            idx += -str;  // =>down
          else
            if ((G[idx+1-str]>=G[idx-str])&&(G[idx+1-str]>=G[idx-1-str]))
              idx += 1-str; // =>down-right
            else
              idx += -1-str;  // =>down-left
        }
      }

    }
}
*/

#define TESTDIR(d) \
if ((!edge[idx+(d)])&&(D[idx+(d)])&&(best_G<G[idx+(d)])&&(idx_pre!=idx+(d))) \
{ \
  best_G = G[idx+(d)]; \
  idx_nex = idx+(d); \
}

staticFunc void ownIPL_trackroutes(
  char *edge,
  int* pnEdges,
  SAnchor* pAncs,
  int nAnc,
  const short* G,
  const char* D,
  int xs,
  int ys)
{
  int cAnc,idx,nEdges=0,idx_pre,idx_nex,best_G;
  SAnchor* pAnc;

    // clear edge image
    memset(edge,0,xs*ys);
    // draw borders to stop edge tracing
    memset(edge,255,xs);
    memset(edge+(ys-1)*xs,255,xs);
    for (idx=1;idx<ys-1;idx++)
      edge[idx*xs] = edge[idx*xs+xs-1] = 255;
    // cycle by anchors
    for (cAnc=0;cAnc<nAnc;cAnc++)
    {
      pAnc = pAncs+cAnc;
      idx_pre = idx = pAnc->idx;
//      if (!D[idx])
  //      continue;   // not an anchor (should never happen indeed)
    //  if (edge[idx])
      //  continue;   // is enumerated already
//      // start 'positive' trace (right or down)
  //    trackrt(edge,xs,idx,G,D, 1);
    //  trackrt(edge,xs,idx,G,D,-1);
      while ((!edge[idx])&&(D[idx]))
      {
        idx_nex = idx;
        best_G = 0;
        TESTDIR(1);
        TESTDIR(-1);
        TESTDIR(xs);
        TESTDIR(-xs);
        TESTDIR(1+xs);
        TESTDIR(1-xs);
        TESTDIR(-1+xs);
        TESTDIR(-1-xs);
        if (!best_G)
          break;
        edge[idx] = 127;
      }
      nEdges++;
    }
    // clear borders
    memset(edge,0,xs);
    memset(edge+(ys-1)*xs,0,xs);
    for (idx=1;idx<ys-1;idx++)
      edge[idx*xs] = edge[idx*xs+xs-1] = 0;
    *pnEdges = nEdges;
}

// edge drawing
//   1) Convolve the image with a separable gaussian filter
//   2) Take the dx and dy the first derivatives using [-1,0,1] and [-1,0,1]^T
//   3) Compute the direction magnitude of gradient
//   4) Find anchors (local maxima)
//   5) Sort anchors
//   6) Trace routes (walk along edges)
MIA_RESULT_CODE IPL_FILT_EdgeDraw_v1(
  unsigned char *edge,        // OUT: canny mask
  const unsigned char *image, // IN:  source image
  short ** ppgx,              // OUT: map of horizontal gradients
  short ** ppgy,              // OUT: map of vertical gradients
  int xs,                     // IN:  image size
  int ys,                     // 
  void* pvBuf,                // IN:  temporary buffer
  int* pnLen,                 // IN/OUT: bytes allocated/used
  const char *nambeg)         // IN:  set this to NULL
{
  char *D,*anchorim;
  short *sm1,*sm2,*gx,*gy,*G;
  int sz,*anchorlist,anchornum,edgenum;
  char nam[FILENAME_MAX];

    if ((xs<=0)||(ys<=0))
      return ERR_GEN_NO_DATA;
    if (pnLen==NULL)
      return ERR_GEN_NULLPOINTER;
if (nambeg)
{
sprintf(nam,"%s_0src.bmp",nambeg);
DBGL_FILE_SaveUint8Image(image,xs,xs,ys,nam);
}
    sz = xs*ys*sizeof(short)+   // smoothed-1
         xs*ys*sizeof(short)+   // smoothed-2
         xs*ys*sizeof(short)+   // gx
         xs*ys*sizeof(short)+   // gy
         xs*ys*sizeof(short)+   // G
         xs*ys              +   // D
         xs*ys              +   // anchors image
         xs*ys*2*sizeof(int);   // anchors list
    if (pvBuf==NULL)
    {
      *pnLen = sz;
      return ERR_OK;
    }
    if (*pnLen<sz)
    {
      *pnLen = sz;
      return ERR_GEN_INSUFFICIENT_BUFFER;
    }
    *pnLen = sz;
    anchorlist = (int*)pvBuf;
    sm1 = (short*)(anchorlist+2*xs*ys);
    sm2 = sm1+xs*ys;
    gx = sm2+xs*ys;
    gy = gx+xs*ys;
    G = gy+xs*ys;
    D = (char*)(G+xs*ys);
    anchorim = D+xs*ys;
    if (ppgx)
      *ppgx = gx;
    if (ppgy)
      *ppgy = gy;
    // perform Haussian smoothing
    IPL_PMP_SmoothX_short(sm1,image,xs,ys);
    IPL_PMP_SmoothY_short(sm2,sm1,xs,ys);
if (nambeg)
{
sprintf(nam,"%s_1smx.bmp",nambeg);
DBGL_FILE_SaveUint16Image((const uint16*)sm1,xs,xs,ys,nam);
sprintf(nam,"%s_1smy.bmp",nambeg);
DBGL_FILE_SaveUint16Image((const uint16*)sm2,xs,xs,ys,nam);
}
    // calculate gradient
    ownIPL_gracal(gx,gy,G,D,sm2,xs,ys,GRATHR);
if (nambeg)
{
sprintf(nam,"%s_2gx.bmp",nambeg);
DBGL_FILE_SaveInt16Image(gx,xs,xs,ys,nam);
sprintf(nam,"%s_2gy.bmp",nambeg);
DBGL_FILE_SaveInt16Image(gy,xs,xs,ys,nam);
sprintf(nam,"%s_2gz.bmp",nambeg);
DBGL_FILE_SaveInt16Image(G,xs,xs,ys,nam);
sprintf(nam,"%s_2r.bmp",nambeg);
DBGL_FILE_SaveInt8Image(D,xs,xs,ys,nam);
}
    // get anchors
    ownIPL_getanchors(anchorlist,&anchornum,anchorim,G,D,xs,ys,1,1);
if (nambeg)
{
sprintf(nam,"%s_3anch.bmp",nambeg);
DBGL_FILE_SaveInt8Image(anchorim,xs,xs,ys,nam);
printf("%d anchors\n",anchornum);
}
    // sort anchors
//    qsort(anchorlist,anchornum,2*sizeof(int),sortanch);
    ownIPL_heapSortAnchsDecreasing((SAnchor*)anchorlist,anchornum);
    // track routes
    ownIPL_trackroutes((char*)edge,&edgenum,(SAnchor*)anchorlist,anchornum,G,D,xs,ys);
if (nambeg)
{
sprintf(nam,"%s_4edges.bmp",nambeg);
DBGL_FILE_SaveInt8Image((const int8*)edge,xs,xs,ys,nam);
printf("%d edges\n",edgenum);
}
    return ERR_OK;
}
