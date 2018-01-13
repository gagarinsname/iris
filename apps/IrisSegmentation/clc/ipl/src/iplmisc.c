/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  21 July 2010                                           */
/*      Revision: 1.0.00                                                 */
/*      Purpose:  miscellaneous functions for IPL                        */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/
#define __FILENUM__ 200 // __FILENUM__TAG200

#include <string.h>
#include <math.h>
#include <memory.h>
#include <malloc.h>
#include "errorcodes.h"
#include "ipl.h"

MIA_RESULT_CODE IPL_GetVersionInfo(
  char* pcVersionString,  // OUT: place reserved for asciiz string
  int pnStringSize)       // IN:  string length
{
  char *str = "Version 3.18.0.142 from 20.09.2014";

    if ((int)strlen(str)>=pnStringSize)
      return ERR_GEN_INSUFFICIENT_BUFFER;
    strcpy(pcVersionString,str);
    return ERR_OK;
}

void IPL_DRAW_Line(
  unsigned char* im,    // image
  unsigned int str,    // stride in bytes
  unsigned int xs,     // image size in pixels
  unsigned int ys,     //
  unsigned int col,    // color
  unsigned int x0,     // line begin
  unsigned int y0,     //
  unsigned int x1,     // line end
  unsigned int y1)     //
{
  double x,y,dx,dy;
  unsigned int i,idx;

    // Sorry for Bresenham. No time to implement it carefully.
    // Make straightforward and stupid.
    dx = (int)(x1)-(int)(x0);
    dy = (int)(y1)-(int)(y0);
    x = sqrt(dx*dx+dy*dy);
    dx /= x;
    dy /= x;
    i = (int)(x+.5);
    x = x0;
    y = y0;
    while (i--)
    {
      if ((x+.5>=0)&&(x+.5<xs)&&(y+.5>=0)&&(y+.5<ys))
      {
        idx = (unsigned int)(y+.5)*str+(unsigned int)(x+.5);
        im[idx  ] = (unsigned char)(col    );
      }
      x += dx;
      y += dy;
    }
}

void ownIPL_FillBorders1B(
  unsigned char* im,
  int xs,
  int ys)
{
  int i;

    // lef & rig
    for (i=1;i<ys-1;i++)
    {
      im[i*xs] = im[i*xs+1];
      im[i*xs+xs-1] = im[i*xs+xs-2];
    }
    // top & bot
    memcpy(im,im+xs,xs);
    memcpy(im+(ys-1)*xs,im+(ys-2)*xs,xs);
}

void ownIPL_FillBorders1B_int(
  int* im,
  int xs,
  int ys)
{
  int i;

    // lef & rig
    for (i=1;i<ys-1;i++)
    {
      im[i*xs] = im[i*xs+1];
      im[i*xs+xs-1] = im[i*xs+xs-2];
    }
    // top & bot
    memcpy(im,im+xs,xs*sizeof(im[0]));
    memcpy(im+(ys-1)*xs,im+(ys-2)*xs,xs*sizeof(im[0]));
}

void ownIPL_FillBordersXB(
  unsigned char* dst,
  int xs,
  int ys,
  int kehw,
  int kehh)
{
  int i,j;

    // sides
    for (i=kehh;i<ys-kehh;i++)
      for (j=0;j<kehw;j++)
      {
        dst[i*xs+j] = dst[i*xs+kehw];
        dst[i*xs+(xs-kehw)+j] = dst[i*xs+(xs-kehw)-1];
      }
    // top & bottom
    for (i=0;i<kehh;i++)
    {
      memcpy(dst+i*xs,dst+kehh*xs,xs);
      memcpy(dst+(ys-kehh+i)*xs,dst+(ys-kehh-1)*xs,xs);
    }
}

void ownIPL_FillBordersXB_int(
  int* dst,
  int xs,
  int ys,
  int kehw,
  int kehh)
{
  int i,j;

    // sides
    for (i=kehh;i<ys-kehh;i++)
      for (j=0;j<kehw;j++)
      {
        dst[i*xs+j] = dst[i*xs+kehw];
        dst[i*xs+(xs-kehw)+j] = dst[i*xs+(xs-kehw)-1];
      }
    // top & bottom
    for (i=0;i<kehh;i++)
    {
      memcpy(dst+i*xs,dst+kehh*xs,xs*sizeof(dst[0]));
      memcpy(dst+(ys-kehh+i)*xs,dst+(ys-kehh-1)*xs,xs*sizeof(dst[0]));
    }
}

MIA_RESULT_CODE IPL_DRAW_MaskFromContour(
  unsigned char** pMask,
  int* pX,
  int* pY,
  int* pW,
  int* pH,
//  const SPolarPointData* pPD_src,
  int xc,
  int yc,
  int Count,
  const int* Radii,
  int isClockWise)
{
  unsigned char* im = NULL;
  int *px=NULL,*py,i,xmin=0,xmax,ymin=0,ymax,w=0,h=0;
  MIA_RESULT_CODE res;

    *pMask = NULL;
    // check arguments
    if (Radii==NULL)
      return ERR_GEN_NULLPOINTER;
    if (Count<=0)
      return ERR_GEN_INVALID_PARAMETER;
    isClockWise = (isClockWise)?(-1):1;
    //== step 1: get coordinates
    // allocate array
    if ((px = (int*)malloc(Count*2*4))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto BRD_Centralize_exit;
    }
    py = px+Count;
    // fill array and determine boundaries
    xmin = ymin = 0x7fffffff;
    xmax = ymax = -0x7fffffff;
    for (i=0;i<(int)(Count);i++)
    {
      px[i] = (int)(xc+Radii[i]*(cos(2*PI*i/Count))+.5);
      py[i] = (int)(yc+Radii[i]*(sin(2*PI*i/Count)*isClockWise)+.5);
      if (xmin>px[i])
        xmin = px[i];
      if (ymin>py[i])
        ymin = py[i];
      if (xmax<px[i])
        xmax = px[i];
      if (ymax<py[i])
        ymax = py[i];
    }
    xmin -= 2;
    xmin &= (~3);
    ymin -= 2;
    ymax += 3;
    xmax += 3;
    xmax = (xmax+3)&(~3);
    for (i=0;i<(int)(Count);i++)
    {
      px[i] -= xmin;
      py[i] -= ymin;
    }
    w = xmax-xmin;
    h = ymax-ymin;
    //== step 2: draw contour
    if ((im = (unsigned char*)malloc(w*h))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto BRD_Centralize_exit;
    }
    memset(im,0,w*h);
    for (i=1;i<(int)(Count);i++)
      IPL_DRAW_Line(im,w,w,h,255,px[i-1],py[i-1],px[i],py[i]);
    IPL_DRAW_Line(im,w,w,h,255,px[Count-1],py[Count-1],px[0],py[0]);
    //== step 3: calculate Mx and My
    res = IPL_CONOB_FloodFill(
      im,              // input image (completely zeroed inside the function)
      w,               // image stride in elements
      w,                // image size
      h,                //
      xc-xmin,                // center to begin fill
      yc-ymin,                //
      0,          // if(TRUE) connect 8 neibours else - 4
      255);             // filling color
    // return
    res = ERR_OK;
    *pMask = im;
    im = NULL;
BRD_Centralize_exit:;
    if (px)
      free(px);
    if (im)
      free(im);
    *pX = xmin;
    *pY = ymin;
    *pW = w;
    *pH = h;
    return res;
}

MIA_RESULT_CODE IPL_DRAW_MaskFromTwoContours(
  unsigned char* mask,  // size = xs*ys
  int xs,
  int ys,
  int isClockwise,
  int isPupilIncl,
  //const SEyeRefinedData* pERD,
  int xc_include,
  int yc_include,
  int Count_include,
  const int* Radii_include,
  int xc_exclude,
  int yc_exclude,
  int Count_exclude,
  const int* Radii_exclude)
{
  MIA_RESULT_CODE res;
  int i,j,x,y,w,h,dx,dy;
  unsigned char *pucMask=NULL;
  //unsigned char* mask = NULL;

    /*if ((mask = (unsigned char*)malloc(xs*ys))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto CalcMaskFromOcclusionRadii_exit;
    }*/
    memset(mask,0,xs*ys);
    // create iris mask
//    if ((res = IPL_DRAW_MaskFromContour(&pucMask,&x,&y,&w,&h,
  //                  &(pERD->iris),isClockwise))!=ERR_OK)
    if ((res = IPL_DRAW_MaskFromContour(&pucMask,&x,&y,&w,&h,
                    xc_include,
                    yc_include,
                    Count_include,
                    Radii_include,
                    isClockwise))!=ERR_OK)
      goto CalcMaskFromOcclusionRadii_exit;
    dx = 0;
    dy = 0;
    if (x<0)
    {
      w += x;
      dx = -x;
      x = 0;
    }
    if (x+w>xs)
    {
      dx = xs-x-w;
      w = xs-x;
    }
    if (y<0)
    {
      dy = -y;
      h += y;
      y = 0;
    }
    if (y+h>ys)
      h = ys-y;
    // copy iris mask to our image-sized mask
    for (i=0;i<h;i++)
      for (j=0;j<w;j++)
        mask[(y+i)*xs+(x+j)] = pucMask[(i+dy)*(w+abs(dx))+(dx>0?dx:0)+j];
    free(pucMask);
    pucMask = NULL;
    if (isPupilIncl)
    {
      // create pupil mask
//      if ((res = IPL_DRAW_MaskFromContour(&pucMask,&x,&y,&w,&h,
  //                &(pERD->pupil),isClockwise))!=ERR_OK)
      if ((res = IPL_DRAW_MaskFromContour(&pucMask,&x,&y,&w,&h,
                  xc_exclude,
                  yc_exclude,
                  Count_exclude,
                  Radii_exclude,
                  isClockwise))!=ERR_OK)
        goto CalcMaskFromOcclusionRadii_exit;
      if (x<0)
      {
        w += x;
        x = 0;
      }
      if (y<0)
      {
        h += y;
        y = 0;
      }
      if (x+w>xs)
        w = xs-x;
      if (y+h>ys)
        h = ys-y;
      // cut away pupil mask
      for (i=0;i<h;i++)
        for (j=0;j<w;j++)
          mask[(y+i)*xs+(x+j)] &= ~(pucMask[i*w+j]);
      free(pucMask);
      pucMask = NULL;
    }
CalcMaskFromOcclusionRadii_exit:;
    if (pucMask)
      free(pucMask);
    //*pmask = mask;
    return res;
}
