/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  11 November 2001                                       */
/*      Modified: 14 February 2002 - adapted IPL_GetOffset               */
/*      Modified: 23 May 2002 - slack code removed (MIA)                 */
/*      Modified: 1 April 2003 - getchesscrosses rewritten (MIA)         */
/*      Modified: 1 July 2003 - Combine function debugged                */
/*      Modified: 12 September 2003 - Pure C code passed to LinesC.c     */
/*      Revision: 1.4.00                                                 */
/*      Purpose:  Searching of straight line crosses on images           */
/*      Authors:                                                         */
/*        Maxim Piskarev                                                 */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#define __FILENUM__ 2 // __FILENUM__TAG2

#include <memory.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "../ownipl.h"
#include "ipl.h"
#include "affine.h"

// sorting by distance
staticFunc int ownIPL_Helper_SortCrossDist(
  const void* e1,
  const void* e2)
{
    return (int)(100.*
      ( (((SDot*)e1)->x)*(((SDot*)e1)->x)+
        (((SDot*)e1)->y)*(((SDot*)e1)->y)-
        (((SDot*)e2)->x)*(((SDot*)e2)->x)-
        (((SDot*)e2)->y)*(((SDot*)e2)->y)) );
//      ( (((SDot*)e1)->x-sg_xc)*(((SDot*)e1)->x-sg_xc)+
  //      (((SDot*)e1)->y-sg_yc)*(((SDot*)e1)->y-sg_yc)-
    //    (((SDot*)e2)->x-sg_xc)*(((SDot*)e2)->x-sg_xc)-
      //  (((SDot*)e2)->y-sg_yc)*(((SDot*)e2)->y-sg_yc)) );
}

// finding best correspondence
staticFunc int ownIPL_ADJ_FindCorr(
  const SDot* base,
  const SDot* scan,
  unsigned int nscan,
  SImageAffineMapping* mapping,
  unsigned int xs,
  unsigned int ys)
{
  double x,y,mindist,dist;
  unsigned int i,j;

    x = base->x;
    y = base->y;
    SImageAffineMapping_Map(mapping,&x,&y);
    if (((int)(x)<0)||((int)(x)>=(int)(xs))||((int)(y)<0)||((int)(y)>=(int)(ys)))
      return -1;  // correspondence point may be outside of image
    mindist = xs*xs+ys*ys;
    j = nscan;
    for (i=0;i<nscan;i++)
      if (base->is13Quadrant==scan[i].is13Quadrant)
      {
        dist = (scan[i].x-x)*(scan[i].x-x)+(scan[i].y-y)*(scan[i].y-y);
        if (dist<mindist)
        {
          j = i;
          mindist = dist;
        }
      }
    if (j==nscan)
      return -2;  // all points are TOO-TOOOOOO FAR
    if (mindist>(xs+ys)*(xs+ys)/3600)
      return -3;  // all points are too far
    return j;
}

// associate crosses in left and right images
MIA_RESULT_CODE IPL_ADJ_CorrespondCrosses2(
  SDot* pLeftCrosses,     // IN/OUT: array of crosses on left image/associated array
  SDot* pRightCrosses,    // IN/OUT: array of crosses on right image/associated array
                          //    items not associated are put in tail
  int* pnumlef, // IN/OUT: number of left crosses/number of black quad 13
  int* pnumrig, // IN/OUT: number of right crosses/number of black quad 24
  int xs,       // IN: image width in pixels
  int ys,       // IN: image height in pixels
  int xl,       // IN: correspondent point coordinates in images
  int yl,       // IN:    (these are used as a reference center)
  int xr,       // IN:
  int yr)       // IN:
{
  int nlef,nrig,cnum,i,j;
  SDot tmpDot;
  SImageAffineMapping ciam_fw,ciam_bw;

    // check arguments
    if ((pLeftCrosses==NULL)||(pRightCrosses==NULL)||(pnumlef==NULL)||(pnumrig==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((xs==0)||(ys==0))
      return ERR_GEN_NO_DATA;
    if ((*pnumlef==0)||(*pnumrig==0))
      return ERR_GEN_INVALID_PARAMETER;
    SImageAffineMapping_SImageAffineMapping(&ciam_fw);
    SImageAffineMapping_SImageAffineMapping(&ciam_bw);
    // sort crosses by decreasing distance from center
//    sg_xc = xl;
  //  sg_yc = yl;
    for (i=0;i<*pnumlef;i++)
    {
      pLeftCrosses[i].x -= xl;
      pLeftCrosses[i].y -= yl;
    }
    qsort(pLeftCrosses,*pnumlef,sizeof(pLeftCrosses[0]),&ownIPL_Helper_SortCrossDist);
    for (i=0;i<*pnumlef;i++)
    {
      pLeftCrosses[i].x += xl;
      pLeftCrosses[i].y += yl;
    }
    //sg_xc = xr;
    //sg_yc = yr;
    for (i=0;i<*pnumrig;i++)
    {
      pRightCrosses[i].x -= xr;
      pRightCrosses[i].y -= yr;
    }
    qsort(pRightCrosses,*pnumrig,sizeof(pRightCrosses[0]),&ownIPL_Helper_SortCrossDist);
    for (i=0;i<*pnumrig;i++)
    {
      pRightCrosses[i].x += xr;
      pRightCrosses[i].y += yr;
    }
    // associate until one of arrays is emptied
    nlef = *pnumlef;
    nrig = *pnumrig;
    cnum = 0;
    if (nlef<nrig)
    { // left image is base one
      SImageAffineMapping_UpdateMapping(&ciam_fw,xl,yl,xr,yr);
      SImageAffineMapping_UpdateMapping(&ciam_bw,xr,yr,xl,yl);
      while (cnum<nlef)
      {
        // find corresponding point
        if ((j = i = ownIPL_ADJ_FindCorr(pLeftCrosses+cnum,pRightCrosses+cnum,
                                            nrig-cnum,&ciam_fw,xs,ys))>=0)
          // reverse pass
          if ((i = ownIPL_ADJ_FindCorr(pRightCrosses+cnum+i,pLeftCrosses+cnum,
                                          nlef-cnum,&ciam_bw,xs,ys))==0)
          { // OK. move correspondent point to appropriate place
            memcpy(&tmpDot,pRightCrosses+cnum+j,sizeof(tmpDot));
            memmove(pRightCrosses+cnum+1,pRightCrosses+cnum,j*sizeof(tmpDot));
            memcpy(pRightCrosses+cnum,&tmpDot,sizeof(tmpDot));
            // update forward and backward mappings
            SImageAffineMapping_UpdateMapping(&ciam_fw,pLeftCrosses [cnum].x,pLeftCrosses [cnum].y,
                                  pRightCrosses[cnum].x,pRightCrosses[cnum].y);
            SImageAffineMapping_UpdateMapping(&ciam_bw,pRightCrosses[cnum].x,pRightCrosses[cnum].y,
                                  pLeftCrosses [cnum].x,pLeftCrosses [cnum].y);
            if ((++cnum)==4)
            {
              SImageAffineMapping_UnupdateMapping(&ciam_fw,xl,yl,xr,yr);
              SImageAffineMapping_UnupdateMapping(&ciam_bw,xr,yr,xl,yl);
            }
            continue;
          }
        // forward or reverse pass failed. exclude the current point
        memmove(pLeftCrosses+cnum,pLeftCrosses+cnum+1,(nlef-cnum-1)*sizeof(tmpDot));
        nlef--;
      }
    }
    else
    { // right image is base one
      SImageAffineMapping_UpdateMapping(&ciam_fw,xr,yr,xl,yl);
      SImageAffineMapping_UpdateMapping(&ciam_bw,xl,yl,xr,yr);
      while (cnum<nrig)
      {
        // find corresponding point
        if ((j = i = ownIPL_ADJ_FindCorr(pRightCrosses+cnum,pLeftCrosses+cnum,
                                            nlef-cnum,&ciam_fw,xs,ys))>=0)
          // reverse pass
          if ((i = ownIPL_ADJ_FindCorr(pLeftCrosses+cnum+i,pRightCrosses+cnum,
                                          nrig-cnum,&ciam_bw,xs,ys))==0)
          { // OK. move correspondent point to appropriate place
            memcpy(&tmpDot,pLeftCrosses+cnum+j,sizeof(tmpDot));
            memmove(pLeftCrosses+cnum+1,pLeftCrosses+cnum,j*sizeof(tmpDot));
            memcpy(pLeftCrosses+cnum,&tmpDot,sizeof(tmpDot));
            SImageAffineMapping_UpdateMapping(&ciam_fw,pRightCrosses[cnum].x,pRightCrosses[cnum].y,
                                  pLeftCrosses [cnum].x,pLeftCrosses [cnum].y);
            SImageAffineMapping_UpdateMapping(&ciam_bw,pLeftCrosses [cnum].x,pLeftCrosses [cnum].y,
                                  pRightCrosses[cnum].x,pRightCrosses[cnum].y);
            if ((++cnum)==4)
            {
              SImageAffineMapping_UnupdateMapping(&ciam_fw,xr,yr,xl,yl);
              SImageAffineMapping_UnupdateMapping(&ciam_bw,xl,yl,xr,yr);
            }
            continue;
          }
        // forward or reverse pass failed. exclude the current point
        memmove(pRightCrosses+cnum,pRightCrosses+cnum+1,(nlef-cnum-1)*sizeof(tmpDot));
        nrig--;
      }
    }
    // count and output number of correspondent quad13 and 24
    for (i=nlef=0;i<(int)(cnum);i++)
      if (pLeftCrosses[i].is13Quadrant)
        nlef++;
    *pnumlef = nlef;
    *pnumrig = cnum-nlef;
    //+++
    // sorting !!!
    //---
    return ERR_OK;
}

// get grosses positions in image of chessboard
// requires four buffers
MIA_RESULT_CODE IPL_ADJ_GetChessCrosses2(
  SDot* pCrosses,           // OUT: crosses found
  int* num,       // IN/OUT: reserved space/number of actually found crosses
  const unsigned char* im,  // IN: source image
  int xs,         // IN: image width in pixels
  int ys,         // IN: image height in pixels
  unsigned char** bufs)     // temporary buffers
{
  enum
  {
    THRESHOLD_1 = 8,
    APERTURE_3 = 6,
    STROKE_LEN = APERTURE_3+1,
  };
  MIA_RESULT_CODE res;
  int objnum_13,objnum_24,i,objnum_out;
  SConObjPrimaryParams *pCOPP_13,*pCOPP_24;
  double x,y;

    // check arguments
    if ((pCrosses==NULL)||(num==NULL)||(im==NULL)||(bufs==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((bufs[0]==NULL)||(bufs[1]==NULL)||(bufs[2]==NULL)||(bufs[3]==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((xs==0)||(ys==0)||(*num==0))
      return ERR_GEN_NO_DATA;
    objnum_out = 0;
    // median filtering with kernel 3x3
//    if ((res = IPL_FilteringMedian(bufs[0],im,GRL_TYPE_BYTE,1,xs,ys,xs,0,0,xs,ys))!=ERR_OK)
    if ((res = IPL_FILT_Median3x3(bufs[0],im,xs,ys))!=ERR_OK)
      return res;
    // allocate for connected objects
    if ((pCOPP_13 = (SConObjPrimaryParams*)malloc(sizeof(SConObjPrimaryParams)*(*num)*2))==NULL)
      return ERR_GEN_NOMEMORY;
    pCOPP_24 = pCOPP_13+(*num);
    // find chessboard crosses with black quadrants 1 and 3
    IPL_FILT_Erode3x3(bufs[1],bufs[0],xs,ys);
    IPL_FILT_Dilate3x3(bufs[2],bufs[0],xs,ys);
    IPL_FILT_ChessCrossesWithBlack13Quadrant(bufs[3],bufs[1],bufs[2],xs,xs,ys,APERTURE_3);
//    memcpy(bufs[5],bufs[3],xs*ys);
    // filter
    for (i=xs*ys;i;i--)
      if (bufs[3][i-1]<THRESHOLD_1)
        bufs[3][i-1] = 0;
      else
        bufs[3][i-1] = (unsigned char)(((int)(bufs[1][i-1])+bufs[2][i-1]+bufs[3][i-1])/3);
    // find connected objects
    objnum_13 = *num;
    if ((res = IPL_CONOB_EnumerateConnectedObjects(
      pCOPP_13,&objnum_13,bufs[3],xs,xs,ys,9,-1,1))!=ERR_OK)
    {
      free(pCOPP_13);
      return res;
    }
    // check, refine, and pack to output
    for (i=0;i<objnum_13;i++)
    {
      x = (double)(pCOPP_13[i].MX)/pCOPP_13[i].M;
      y = (double)(pCOPP_13[i].MY)/pCOPP_13[i].M;
      if ((x>STROKE_LEN)&&(x<xs-STROKE_LEN-1)&&(y>STROKE_LEN)&&(y<ys-STROKE_LEN-1))
//        if (ownIPL_ADJ_RefineCross(&x,&y,im,xs)<(50<<16))
//        if (ownIPL_ADJ_CheckCross(&x,&y,im,xs,1)==ERR_OK)
        { // ok. pack it
          pCrosses[objnum_out].x = x;
          pCrosses[objnum_out].y = y;
          pCrosses[objnum_out].is13Quadrant = 1;
          pCrosses[objnum_out].isDefined = 1;
          if ((int)(++objnum_out)>=(*num))
          {
            free(pCOPP_13);
            return ERR_GEN_INSUFFICIENT_BUFFER;
          }
        }
    }
    // find chessboard crosses with black quadrants 2 and 4
//    IPL_FilteringMorphFixed(bufs[1],bufs[0],GRL_TYPE_BYTE,GRL_MORPH_ERODE,GRL_MASK_RECTANGLE,xs,ys,xs,0,0,xs,ys);
//    IPL_FilteringMorphFixed(bufs[2],bufs[0],GRL_TYPE_BYTE,GRL_MORPH_DILATE,GRL_MASK_RECTANGLE,xs,ys,xs,0,0,xs,ys);
    IPL_FILT_ChessCrossesWithBlack24Quadrant(bufs[3],bufs[1],bufs[2],xs,xs,ys,APERTURE_3);
//    memcpy(bufs[5],bufs[3],xs*ys);
    // filter
    for (i=xs*ys;i;i--)
      if (bufs[3][i-1]<THRESHOLD_1)
        bufs[3][i-1] = 0;
      else
        bufs[3][i-1] = (unsigned char)(((int)(bufs[1][i-1])+bufs[2][i-1]+bufs[3][i-1])/3);
    // find connected objects
    objnum_24 = *num;
    if ((res = IPL_CONOB_EnumerateConnectedObjects(
      pCOPP_24,&objnum_24,bufs[3],xs,xs,ys,9,-1,1))!=ERR_OK)
    {
      free(pCOPP_13);
      return res;
    }
    // check, refine, and pack to output
    for (i=0;i<objnum_24;i++)
    {
      x = (double)(pCOPP_24[i].MX)/pCOPP_24[i].M;
      y = (double)(pCOPP_24[i].MY)/pCOPP_24[i].M;
      if ((x>STROKE_LEN)&&(x<xs-STROKE_LEN-1)&&(y>STROKE_LEN)&&(y<ys-STROKE_LEN-1))
//        if (ownIPL_ADJ_RefineCross(&x,&y,im,xs)<(50<<16))
//        if (ownIPL_ADJ_CheckCross(&x,&y,im,xs,0)==ERR_OK)
        { // ok. pack it
          pCrosses[objnum_out].x = x;
          pCrosses[objnum_out].y = y;
          pCrosses[objnum_out].is13Quadrant = 0;
          pCrosses[objnum_out].isDefined = 1;
          if ((int)(++objnum_out)>=(*num))
          {
            free(pCOPP_13);
            return ERR_GEN_INSUFFICIENT_BUFFER;
          }
        }
    }
    *num = objnum_out;
    free(pCOPP_13);
    return ERR_OK;
}

// subpixel refinement of points
staticFunc MIA_RESULT_CODE ownIPL_ADJ_RefineOneCross(
  SDot* pCross,             // OUT: cross
  const unsigned char* im,  // IN: source image
  int str,        // IN: image stride
  int w,          // IN: window size in pixels
  int h)          //     caller guarantees window does not span out of image
{
  double a[4],b[2],d;
  int gx,gy,i,j;

    // initialise linear system
    memset(a,0,sizeof(a));
    memset(b,0,sizeof(b));
    // sum up in region
    for (i=h-2;i;i--)
      for (j=w-2;j;j--)
      {
        gx = (int)(im[i*str+j+1])-(int)(im[i*str+j-1]);
        gy = (int)(im[(i+1)*str+j])-(int)(im[(i-1)*str+j]);
        a[0] += gx*gx;
        a[1] += gx*gy;
        a[2] += gx*gy;
        a[3] += gy*gy;
        b[0] += j*gx*gx+i*gx*gy;
        b[1] += j*gx*gy+i*gy*gy;
      }
    // calculate coefficients
    if ((d = a[0]*a[3]-a[1]*a[2])==0.)
      return ERR_IPL_UNEXPRESSIVE_DATA;
    pCross->x = (b[0]*a[3]-b[1]*a[1])/d;
    pCross->y = (b[1]*a[0]-b[0]*a[2])/d;
    return ERR_OK;
}

// subpixel refinement of points
MIA_RESULT_CODE IPL_ADJ_RefineChessCrosses(
  SDot* pCrosses,           // OUT: crosses
  int num,        // IN: number of crosses
  const unsigned char* im,  // IN: source image
  int xs,         // IN: image width in pixels
  int ys,         // IN: image height in pixels
  int hws,        // IN: recommened horizontal half-size of window
  int vws)        // IN: recommened vertical half-size of window
{                           //     if zero, functions guesses by itself
  int i,x,y;
  SDot Cross;

    // check arguments
    if ((pCrosses==NULL)||(im==NULL))
      return ERR_GEN_NULLPOINTER;
    // determine whs
    if (!hws)
      // 0.034 is a magic number. You not shaman - you not touch!
      hws = vws = (int)(sqrt(.034*xs*ys/num)+.5);
    else
      // set outside, check reasonability
      if ((xs<hws*2+1)||(ys<vws*2+1))
        return ERR_GEN_NO_DATA;
    // cycle by crosses
    for (i=0;i<num;i++)
    {
      // get integer position
      x = (int)(pCrosses[i].x+.5);
      y = (int)(pCrosses[i].y+.5);
      // check for outspanning
      if (x<hws)
        x = hws;
      if (y<vws)
        y = vws;
      if (x+hws>=xs)
        x = xs-hws-1;
      if (y+vws>=ys)
        y = ys-vws-1;
      // refine
      if (ownIPL_ADJ_RefineOneCross(&Cross,
            im+(y-vws)*xs+(x-hws),xs,
            2*hws+1,2*vws+1)==ERR_OK)
      { // correct position
        pCrosses[i].x = x+Cross.x-hws;
        pCrosses[i].y = y+Cross.y-vws;
      }
    }
    return ERR_OK;
}

typedef struct
{
  SDot Point;       // x,y, and quadrant flag
  int Neighbors[4]; // neighbour number in the array
  int posx;         // position in mesh relative to beginner point
  int posy;         //   beginner point has (0,0)
} SDotEx;

staticFunc int ownIPL_ADJ_SortSDotEx(
  const void* e1,
  const void* e2)
{
    if (((SDotEx*)e1)->posy!=((SDotEx*)e2)->posy)
      return ((SDotEx*)e1)->posy-((SDotEx*)e2)->posy;
    return ((SDotEx*)e1)->posx-((SDotEx*)e2)->posx;
}

// ordering points in a rectangular mesh
staticFunc MIA_RESULT_CODE IPL_ADJ_OrderPointsInMesh(
  SDot** outcrosses,    // OUT: outputted points in a raster (rectangle mesh) form
  const SDot* crosses,  // IN:  input points array
  int* num,   // IN:  number of points, OUT: number of approved points
  int* nx,    // OUT: horizontal size of mesh
  int* ny,    // OUT: vertical size of mech
  int xs,     // IN:  image width
  int ys,     // IN:  image height
  int* pMeshSize) // OUT: average size of grid cell in pixels
{
  MIA_RESULT_CODE retval = ERR_GEN_NOMEMORY;
  enum
  {
    NOT_BURNED = 0x80000000,    // signs dot was not touched by forest fire
    ENOUGH_DOTS = 3,
    ANGLE_JITTER = 15,
    GRIDJITTER = 4, // jitter of the grid
  };
  SDot* out = NULL;
  SDotEx *DotList = NULL;
  int *histogram = NULL,dummy_nx,dummy_ny;
  int nq13,nq24,diag,i,j,k,MeshSize,MeshTilt;
  int mx,Mx,my,My,mx_,Mx_,my_,My_;

    // check arguments
    if ((outcrosses==NULL)||(crosses==NULL)||(num==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((xs==0)||(ys==0)||(*num==0))
      return ERR_GEN_NO_DATA;
    if (nx==NULL)
      nx = &dummy_nx;
    if (ny==NULL)
      ny = &dummy_ny;
    // calculate number of quadrant13 and quadrant24 points
    nq13 = nq24 = 0;
    for (i=0;i<*num;i++)
      if (crosses[i].is13Quadrant)
        nq13++;
      else
        nq24++;
    // test if enough points of both types
    if ((nq13<(*num/4))||(nq24<(*num/4)))
      return ERR_IPL_SPARSE_POINTS;
    // calculte image diagonal (biggest distance possible
    diag = (int)(sqrt(xs*xs+ys*ys)+1);
    // allocate histogram array
    i = (diag>90)?diag:90;
    i = (i>*num)?i:*num;
    if ((histogram = (int*)malloc(sizeof(int)*i))==NULL)
      goto opim_exit;
    // clear histogram
    memset(histogram,0,sizeof(histogram[0])*diag);
    // estimation of mesh unit size
    // fill histogram of all inter-point distances (for points with different quadrant)
    for (i=0;i<*num;i++)
      for (j=0;j<i;j++)
        if (crosses[j].is13Quadrant!=crosses[i].is13Quadrant)
        {
          k = (int)
            sqrt((crosses[i].x-crosses[j].x)*(crosses[i].x-crosses[j].x)+
                 (crosses[i].y-crosses[j].y)*(crosses[i].y-crosses[j].y));
          histogram[k]++;
          // blur histogram to eliminate jitter
          if (k)
            histogram[k-1]++;
          if (k+1<diag)
            histogram[k+1]++;
        }
    // locate first pronounced maxima
    for (i=0;i<diag;i++)
      if (histogram[i]>(*num)/2)
        break;
    if (i==diag)
    { // all histogram searched, sufficient maxima not discovered
      retval = ERR_GEN_INVALID_PARAMETER;
      goto opim_exit;
    }
    // refine maxima (we found first outspan, peak maybe higher)
    while ((i<diag-3)&&((histogram[i+1]>histogram[i])||(histogram[i+2]>histogram[i])))
      // each step looks forward up to two points
      if (histogram[i+1]>histogram[i])
        i += 1;
      else
        i += 2;
    if (i>=diag-4)
    { // very strange. optimised to a boundary
      retval = ERR_GEN_INTERNAL;
      goto opim_exit;
    }
    // remember this as the size of mesh
    MeshSize = i;
    // estimation of tilt angle in degrees
    // clear the histogram
    memset(histogram,0,sizeof(histogram[0])*90);
    // fill the angular histogram (for points with different quadrant)
    for (i=0;i<*num;i++)
      for(j=0;j<i;j++)
        if ((crosses[j].is13Quadrant!=crosses[i].is13Quadrant)&&
            (abs((int)(GRIDJITTER*(sqrt((crosses[i].x-crosses[j].x)*
                                       (crosses[i].x-crosses[j].x)+
                                       (crosses[i].y-crosses[j].y)*
                                       (crosses[i].y-crosses[j].y))-MeshSize)))<(int)(MeshSize)))
        {
          k = (int)((PI+atan2((crosses[i].y-crosses[j].y),
                              (crosses[i].x-crosses[j].x)))*180/PI);
          k %= 90;
          // blur to prevent jittering
          histogram[k]++;
          histogram[(k+1 )%90]++;
          histogram[(k+89)%90]++;
        }
    // find dominating direction
    for (i=0;i<90;i++)
      if (histogram[i]>*num)
        break;
    if (i==90)
    { // all histogram searched, sufficient maxima not discovered
      retval = ERR_GEN_BAD_ALIGNMENT;
      goto opim_exit;
    }
    // refine maxima (we found first outspan, peak maybe higher)
    j = i;
    while ((histogram[(i+1)%90]>histogram[i])||(histogram[(i+2)%90]>histogram[i]))
      // each step looks forward up to two points
      if (histogram[(i+1)%90]>histogram[i])
        i = (i+1)%90;
      else
        i = (i+2)%90;
    // since circled, we can go reverce
    while ((histogram[(j+89)%90]>histogram[j])||(histogram[(j+88)%90]>histogram[j]))
      // each step looks forward up to two points
      if (histogram[(j+89)%90]>histogram[j])
        j = (j+89)%90;
      else
        j = (j+88)%90;
    MeshTilt = (histogram[i]>histogram[j])?i:j;
    // if upper octant, drop below zero
    if (MeshTilt>45)
      MeshTilt += 270;
    // enumerate neighbors of each point
    if ((DotList = (SDotEx*)malloc(sizeof(SDotEx)*(*num)))==NULL)
      goto opim_exit;
    // fill list of dots with neighbors
    for (i=0;i<*num;i++)
    {
      // copy data from source point
      memcpy(&DotList[i].Point,crosses+i,sizeof(crosses[0]));
      // no neighbors traced yet
      DotList[i].Neighbors[0] = -1;
      DotList[i].Neighbors[1] = -1;
      DotList[i].Neighbors[2] = -1;
      DotList[i].Neighbors[3] = -1;
      // at the same time, initialise for forest fire algorithm
      DotList[i].posx = NOT_BURNED;
      // try each point for neighbor (with other quadrant)
      for (j=0;j<*num;j++)
        if ((crosses[j].is13Quadrant!=crosses[i].is13Quadrant)&&
            (abs((int)(GRIDJITTER*(sqrt((crosses[i].x-crosses[j].x)*
                                       (crosses[i].x-crosses[j].x)+
                                       (crosses[i].y-crosses[j].y)*
                                       (crosses[i].y-crosses[j].y))-MeshSize)))<(int)(MeshSize)))
        {
          k = (int)((PI+atan2((crosses[i].y-crosses[j].y),
                              (crosses[i].x-crosses[j].x)))*180/PI);
          // east neighbor?
          if (((k-MeshTilt+540)%360>180-ANGLE_JITTER)&&((k-MeshTilt+540)%360<180+ANGLE_JITTER))
            DotList[i].Neighbors[0] = j;
          // north neighbor?
          if (((k-MeshTilt+630)%360>180-ANGLE_JITTER)&&((k-MeshTilt+630)%360<180+ANGLE_JITTER))
            DotList[i].Neighbors[1] = j;
          // west neighbor?
          if (((k-MeshTilt+720)%360>180-ANGLE_JITTER)&&((k-MeshTilt+720)%360<180+ANGLE_JITTER))
            DotList[i].Neighbors[2] = j;
          // south neighbor?
          if (((k-MeshTilt+810)%360>180-ANGLE_JITTER)&&((k-MeshTilt+810)%360<180+ANGLE_JITTER))
            DotList[i].Neighbors[3] = j;
        }
    }
    // look for point with four (or at least three) neighbours
    for (i=0;i<*num;i++)
      if ((DotList[i].Neighbors[0]!=-1)&&
          (DotList[i].Neighbors[1]!=-1)&&
          (DotList[i].Neighbors[2]!=-1)&&
          (DotList[i].Neighbors[3]!=-1))
        break;
    if (i==*num)
    { // sorry, no points have four neighbours. maybe, three?
      for (i=0;i<*num;i++)
        if (((DotList[i].Neighbors[0]!=-1)?1:0)+
            ((DotList[i].Neighbors[1]!=-1)?1:0)+
            ((DotList[i].Neighbors[2]!=-1)?1:0)+
            ((DotList[i].Neighbors[3]!=-1)?1:0) == 3)
          break;
      if (i==*num)
      { // very sparse mesh. unacceptable
        retval = ERR_IPL_SPARSE_POINTS;
        goto opim_exit;
      }
    }
    // starting point found. use forest fire now
    DotList[i].posx = 0;
    DotList[i].posy = 0;
    // histogram now serves as LIFO index buffer
    histogram[j = 0] = i;
    j++;
    while (j)
    {
      i = histogram[--j];
      // east neighbor
      if (DotList[k = DotList[i].Neighbors[0]].posx==NOT_BURNED)
      {
        DotList[k].posx = DotList[i].posx+1;
        DotList[k].posy = DotList[i].posy;
        histogram[j++] = k;
      }
      // north neighbor
      if (DotList[k = DotList[i].Neighbors[1]].posx==NOT_BURNED)
      {
        DotList[k].posx = DotList[i].posx;
        DotList[k].posy = DotList[i].posy-1;
        histogram[j++] = k;
      }
      // west neighbor
      if (DotList[k = DotList[i].Neighbors[2]].posx==NOT_BURNED)
      {
        DotList[k].posx = DotList[i].posx-1;
        DotList[k].posy = DotList[i].posy;
        histogram[j++] = k;
      }
      // south neighbor
      if (DotList[k = DotList[i].Neighbors[3]].posx==NOT_BURNED)
      {
        DotList[k].posx = DotList[i].posx;
        DotList[k].posy = DotList[i].posy+1;
        histogram[j++] = k;
      }
    }
    // search for minimal and maximal x and y
    for (i=mx=my=Mx=My=0;i<*num;i++)
      if (DotList[i].posx!=NOT_BURNED)
      {
        if (mx>DotList[i].posx)
          mx = DotList[i].posx;
        if (my>DotList[i].posy)
          my = DotList[i].posy;
        if (Mx<DotList[i].posx)
          Mx = DotList[i].posx;
        if (My<DotList[i].posy)
          My = DotList[i].posy;
      }
    // make a histogram of x's (shifted by min(x))
    memset(histogram,0,*num*sizeof(histogram[0]));
    for (i=0;i<*num;i++)
      if (DotList[i].posx!=NOT_BURNED)
        histogram[DotList[i].posx-mx]++;
    // see first position with not less than ENOUGH_DOTS points of mesh
    for (i=0;i<Mx-mx;i++)
      if (histogram[i]>=ENOUGH_DOTS)
        break;
    if (i==(Mx-mx))
    { // no vertical line of mesh contains ENOUGH_DOTS points. bad (unexpressive) mesh
      retval = ERR_IPL_UNEXPRESSIVE_DATA;
      goto opim_exit;
    }
    // set zero level of x to position with not less than ENOUGH_DOTS points of mesh
    mx_ = mx+i;
    // see last position with not less than ENOUGH_DOTS points of mesh
    for (i=Mx+1-mx;i;i--)
      if (histogram[i-1]>=ENOUGH_DOTS)
        break;
    // set max level of x to position with not less than ENOUGH_DOTS points of mesh
    Mx_ = mx+i;
    // make a histogram of x's (shifted by min(x))
    memset(histogram,0,*num*sizeof(histogram[0]));
    for (i=0;i<*num;i++)
      if (DotList[i].posx!=NOT_BURNED)
        histogram[DotList[i].posy-my]++;
    // see first position with not less than ENOUGH_DOTS points of mesh
    for (i=0;i<(My-my);i++)
      if (histogram[i]>=ENOUGH_DOTS)
        break;
    if (i==(My-my))
    { // no horizontal line of mesh contains ENOUGH_DOTS points. bad (unexpressive) mesh
      retval = ERR_IPL_UNEXPRESSIVE_DATA;
      goto opim_exit;
    }
    // set zero level of y to position with not less than ENOUGH_DOTS points of mesh
    my_ = my+i;
    // see last position with not less than ENOUGH_DOTS points of mesh
    for (i=My+1-my;i;i--)
      if (histogram[i-1]>=ENOUGH_DOTS)
        break;
    // set max level of y to position with not less than ENOUGH_DOTS points of mesh
    My_ = my+i;
    // allocate array of dots of required dimension
    if ((out = (SDot*)malloc(sizeof(SDot)*(Mx_-mx_)*(My_-my_)))==NULL)
      goto opim_exit;
    memset(out,0,sizeof(SDot)*(Mx_-mx_)*(My_-my_));
    // sort points. this makes rasterization
    qsort(DotList,*num,sizeof(DotList[0]),&ownIPL_ADJ_SortSDotEx);
    // output points
    for (j=i=0;i<*num;i++)
      if ((DotList[i].posx!=NOT_BURNED)&&
          (DotList[i].posx>=mx_)&&(DotList[i].posx<Mx_)&&
          (DotList[i].posy>=my_)&&(DotList[i].posy<My_) )
      {
        j++;
        k = (DotList[i].posy-my_)*(Mx_-mx_)+(DotList[i].posx-mx_);
        out[k].x = DotList[i].Point.x ;
        out[k].y = DotList[i].Point.y;
        out[k].is13Quadrant = DotList[i].Point.is13Quadrant;
        out[k].isDefined = 1;
      }
    *outcrosses = out;
    // output dimensions
    *num = j;         // OUT: number of approved points
    *nx = Mx_-mx_;    // OUT: horizontal size of mesh
    *ny = My_-my_;    // OUT: vertical size of mech
    retval = ERR_OK;
    if(pMeshSize)
      *pMeshSize = MeshSize;
opim_exit:
    if (histogram)
      free(histogram);
    if (DotList)
      free(DotList);
    return retval;
}

typedef struct
{
  double quality;     // deviation from straightness (sum of point distances to the line)
  double uncertainty; // uncertainty of quality detection
  double radius;      // coordinates in radial system: length of perpendicular
  double slope;       // coordinates in radial system: angle of perpendicular
  int len;  // how many points in line
} SMeshLineParameters;

// calculate the value showing deviation of line from ideally straight
staticFunc void IPL_ADJ_EvaluateLineQuality(
  SMeshLineParameters* mlp, // OUT: calculated parameters
  const SDot* mesh,         // IN:  start of line in mesh
  int str,        // IN:  stride between line elements
  int n,          // IN:  length of line
  int *pidx)      // TMP: index buffer of n elements
{
  double Mx,My,Mxx,Mxy,Myy,rx,ry,D,L1;
  int i,M;
  const SDot* ptr;

    // calculate moments of orders 0,1,2
    for (Mx = My = Mxx = Mxy = Myy = 0.,M = 0,ptr = mesh,i = n;i;i--,ptr += str)
      if (ptr->isDefined)
      {
        pidx[M++] = (int)(mesh-ptr);
        Mx += ptr->x;
        My += ptr->y;
        Mxx += ptr->x*ptr->x;
        Mxy += ptr->x*ptr->y;
        Myy += ptr->y*ptr->y;
      }
    // save length of line
    mlp->len = M;
    if (mlp->len<=1)
    { // single point or no points. No line slope, no quality.
      mlp->quality = FLT_MAX;
      return;
    }
    if (mlp->len==2)
    { // two points. Ideal quality, definite slope.
      mlp->quality = 0.;
      mlp->uncertainty = 0.;
// radius-vector to closest line point for line of two points is calculated as:
//
//     1      s * d
// r = - (s-d -----)
//     2      d * d
//
// where s is sum of point's radius-vectors, d is difference of points radius-vectors
      // let s=(Mx,My), d=(Mxx,Myy)
      Mx = mesh[pidx[0]].x+mesh[pidx[1]].x;
      My = mesh[pidx[0]].y+mesh[pidx[1]].y;
      Mxx = mesh[pidx[0]].x-mesh[pidx[1]].x;
      Myy = mesh[pidx[0]].y-mesh[pidx[1]].y;
      // let Mxy be (s*d)/(d*d)
      Mxy = (Mx*Mxx+My*Myy)/(Mxx*Mxx+Myy*Myy);
      rx = Mx-Mxx*Mxy;
      ry = My-Myy*Mxy;
      mlp->radius = sqrt(rx*rx+ry*ry)/2.;
      mlp->slope = atan2(ry,rx);
      return;
    }
    // finished with dummies. real work goes here.
    // normalise 1st order moments
    Mx /= M;
    My /= M;
    // normalise 2nd order moments
    Mxx = Mxx/M-Mx*Mx;
    Mxy = Mxy/M-Mx*My;
    Myy = Myy/M-My*My;
    // solve eigen equation for 2x2 matrix (Mxx, Mxy, Mxy, Myy)
    //  and get smaller eigenvalue and corresponding eigenvector
    // however, better to get bigger eigenvector (more precision, less exclusions)
    D = sqrt((Mxx-Myy)*(Mxx-Myy)+4*Mxy*Mxy);
    L1 = (Mxx+Myy+D)/2.;
    // norma
    D = sqrt((Mxx-L1)*(Mxx-L1)+Mxy*Mxy);
    // (rx,ry): second eigen
    rx = (Mxx-L1)/D;
    ry = Mxy/D;
    // sum up projections
    for (Mxy = 0.,ptr = mesh,i = n;i;i--,ptr += str)
      if (ptr->isDefined)
      {
        D = rx*(ptr->x-Mx)+ry*(ptr->y-My);
        Mxy += D*D;
      }
    // output
    mlp->quality = Mxy/M;
    mlp->uncertainty = 1./M;
    D = Mx*rx+My*ry;
    mlp->radius = fabs(D);
    mlp->slope = atan2(ry,rx);
    if (D<0.)
    {
      if (mlp->slope>0.)
        mlp->slope -= PI;
      else
        mlp->slope += PI;
    }
}

// fills gaps in mesh with interpolation from existing nodes
staticFunc int ownIPL_ADJ_InterpolateWithGivenNeighborCount(
  SDot* dst,          // OUT: filled mesh
  int nx,   // IN:  mesh dimensions
  int ny,
  int c)    // IN:  required neighbor count
{
  int filled_nodes,i,j,neib_c;
  double xm1,x0,xp1,ym1,y0,yp1;
  int nxm1,nx0,nxp1,nym1,ny0,nyp1;

    filled_nodes = 0;
    // search for empty nodes
    for (i=0;i<ny;i++)
      for (j=0;j<nx;j++)
        if (!dst[i*nx+j].isDefined)
        { // found empty node.
          // count existing neighbors
          xm1 = ym1 = x0 = y0 = xp1 = yp1 = 0.;
          nxm1 = nx0 = nxp1 = nym1 = ny0 = nyp1 = neib_c = 0;
          // upper-left neighbor
          if ((i     )&&(j     )&&(dst[(i-1)*nx+(j-1)].isDefined))
          {
            nxm1++;
            nym1++;
            xm1 += dst[(i-1)*nx+(j-1)].x;
            ym1 += dst[(i-1)*nx+(j-1)].y;
            neib_c++;
          }
          // upper neighbor
          if ((i     )&&          (dst[(i-1)*nx+(j  )].isDefined))
          {
            nx0++;
            nym1++;
            x0 += dst[(i-1)*nx+(j  )].x;
            ym1 += dst[(i-1)*nx+(j  )].y;
            neib_c++;
          }
          // upper-right neighbor
          if ((i     )&&(j+1<nx)&&(dst[(i-1)*nx+(j+1)].isDefined))
          {
            nxp1++;
            nym1++;
            xp1 += dst[(i-1)*nx+(j+1)].x;
            ym1 += dst[(i-1)*nx+(j+1)].y;
            neib_c++;
          }
          // lower-left neighbor
          if ((i+1<ny)&&(j     )&&(dst[(i+1)*nx+(j-1)].isDefined))
          {
            nxm1++;
            nyp1++;
            xm1 += dst[(i+1)*nx+(j-1)].x;
            yp1 += dst[(i+1)*nx+(j-1)].y;
            neib_c++;
          }
          // lower neighbor
          if ((i+1<ny)&&          (dst[(i+1)*nx+(j  )].isDefined))
          {
            nx0++;
            nyp1++;
            x0 += dst[(i+1)*nx+(j  )].x;
            yp1 += dst[(i+1)*nx+(j  )].y;
            neib_c++;
          }
          // lower-right neighbor
          if ((i+1<ny)&&(j+1<nx)&&(dst[(i+1)*nx+(j+1)].isDefined))
          {
            nxp1++;
            nyp1++;
            xp1 += dst[(i+1)*nx+(j+1)].x;
            yp1 += dst[(i+1)*nx+(j+1)].y;
            neib_c++;
          }
          // left neighbor
          if (          (j     )&&(dst[(i  )*nx+(j-1)].isDefined))
          {
            nxm1++;
            ny0++;
            xm1 += dst[(i  )*nx+(j-1)].x;
            y0 += dst[(i  )*nx+(j-1)].y;
            neib_c++;
          }
          // right neighbor
          if (          (j+1<nx)&&(dst[(i  )*nx+(j+1)].isDefined))
          {
            nxp1++;
            ny0++;
            xp1 += dst[(i  )*nx+(j+1)].x;
            y0 += dst[(i  )*nx+(j+1)].y;
            neib_c++;
          }
          // decide if we can fill it
          if ((neib_c>=c) &&
              ( (nx0) || ((nxm1)&&(nxp1)) ) &&
              ( (ny0) || ((nym1)&&(nyp1)) ) )
          { // yes we should fill it now
            if ((!nxm1)||(!nxp1))
              dst[i*nx+j].x = x0/nx0;
            else
              if (!nx0)
                dst[i*nx+j].x = (xp1/nxp1+xm1/nxm1)/2.;
              else
                dst[i*nx+j].x = (x0/nx0+(xp1/nxp1+xm1/nxm1)/2.)/2.;
            if ((!nym1)||(!nyp1))
              dst[i*nx+j].y = y0/ny0;
            else
              if (!ny0)
                dst[i*nx+j].y = (yp1/nyp1+ym1/nym1)/2.;
              else
                dst[i*nx+j].y = (y0/ny0+(yp1/nyp1+ym1/nym1)/2.)/2.;
            // denote the fact
            filled_nodes++;
            dst[i*nx+j].isDefined = 1;
          }
        }
    return filled_nodes;
}

// fills gaps in mesh with interpolation from existing nodes
staticFunc MIA_RESULT_CODE IPL_ADJ_InterpolateCrossesByMesh(
  SDot* dst,          // OUT: filled mesh
  const SDot* src,    // IN:  source mesh
  int nx,   // IN:  mesh dimensions
  int ny)
{
  int i,j,emptynodes;

    // check arguments
    if ((dst==NULL)||(src==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((nx==0)||(ny==0))
      return ERR_GEN_NO_DATA;
    // copy src to dst
    memcpy(dst,src,nx*ny*sizeof(dst[0]));
    // count empty nodes
    for (emptynodes = 0,i = nx*ny;i;i--)
      if (!dst[i-1].isDefined)
        emptynodes++;
    // try filling with 7,6,5,4,3 neighbors in order unless all filled
    while (emptynodes)
    {
      emptynodes -= (j = ownIPL_ADJ_InterpolateWithGivenNeighborCount(dst,nx,ny,7));
      if (j)
        continue;
      emptynodes -= (j = ownIPL_ADJ_InterpolateWithGivenNeighborCount(dst,nx,ny,6));
      if (j)
        continue;
      emptynodes -= (j = ownIPL_ADJ_InterpolateWithGivenNeighborCount(dst,nx,ny,5));
      if (j)
        continue;
      emptynodes -= (j = ownIPL_ADJ_InterpolateWithGivenNeighborCount(dst,nx,ny,4));
      if (j)
        continue;
      emptynodes -= (j = ownIPL_ADJ_InterpolateWithGivenNeighborCount(dst,nx,ny,3));
      if (!j)
        // some points are unreachable even with 3 neighbors.
        // this means whole row or column (or more) is empty.
        // too bad for us.
        return ERR_CLB_TOO_SPARSE_MESH;
    }
    return ERR_OK;
}

// simply get a mean color of small neighborhood of a point
staticFunc int ownIPL_ADJ_GetLocalColor(
  const unsigned char* im,    // IN:  source image
  int xs,           // IN:  image size
  int ys,           // IN:
  int x,                     // IN:  position
  int y)                     // IN:
{
  enum
  {
    NB_SIZE = 3,
  };
  int i,j,S;

    if ((x-NB_SIZE<0)||(y-NB_SIZE<0)||(x+NB_SIZE>=(int)(xs))||(y+NB_SIZE>=(int)(ys)))
      return -1;
    for (S=0,i=-NB_SIZE;i<=NB_SIZE;i++)
      for (j=-NB_SIZE;j<=NB_SIZE;j++)
        S += im[(y+i)*xs+(x+j)];
    return S/((NB_SIZE*2+1)*(NB_SIZE*2+1));
}

// find a white circle in a chessboard
staticFunc MIA_RESULT_CODE IPL_ADJ_FindWhiteCircleMarker(
  SDot* Marker,             // OUT: marker position
  const SDot* crosses,      // IN:  crosses array
  int nx,         // IN:  crosses array size
  int ny,         //
  const unsigned char* im,  // IN:  image
  int xs,         // IN:  image size
  int ys)         //
{
  int i,j,sb,nb,sw,nw;
  int x,y,cc,xdx,xdy,ydx,ydy;
  SDot* pcr;
  MIA_RESULT_CODE res;

    // check arguments
    if ((crosses==NULL)||(im==NULL)||(Marker==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((xs==0)||(ys==0)||(nx==0)||(ny==0))
      return ERR_GEN_NO_DATA;
    // allocate memory
    if ((pcr = (SDot*)malloc(sizeof(SDot)*nx*ny))==NULL)
      return ERR_GEN_NOMEMORY;
    // fill gaps
    if ((res = IPL_ADJ_InterpolateCrossesByMesh(pcr,crosses,nx,ny))!=ERR_OK)
    {
      free(pcr);
      return res;
    }
    // search for marker
    for (i=0;i<ny-1;i++)
      for (j=0;j<nx-1;j++)
        if (!pcr[i*nx+j].is13Quadrant)
        {
          x = (int)((pcr[(i  )*nx+(j  )].x+pcr[(i+1)*nx+(j  )].x+
                    pcr[(i  )*nx+(j+1)].x+pcr[(i+1)*nx+(j+1)].x)/4.+.5);
          y = (int)((pcr[(i  )*nx+(j  )].y+pcr[(i+1)*nx+(j  )].y+
                    pcr[(i  )*nx+(j+1)].y+pcr[(i+1)*nx+(j+1)].y)/4.+.5);
          if ((x>=0)&&(y>=0)&&(x<(int)(xs))&&(y<(int)(ys)))
          {
            // get what is supposed to be white and what is supposed to be black
            // shift of x coordinate for one horizontal mesh step
            xdx = (int)((pcr[(i+1)*nx+(j+1)].x-pcr[(i+1)*nx+(j  )].x+
                        pcr[(i  )*nx+(j+1)].x-pcr[(i  )*nx+(j  )].x)/2.+.5);
            // shift of y coordinate for one horizontal mesh step
            ydx = (int)((pcr[(i+1)*nx+(j+1)].y-pcr[(i+1)*nx+(j  )].y+
                        pcr[(i  )*nx+(j+1)].y-pcr[(i  )*nx+(j  )].y)/2.+.5);
            // shift of x coordinate for one vertical mesh step
            xdy = (int)((pcr[(i+1)*nx+(j+1)].x-pcr[(i  )*nx+(j+1)].x+
                        pcr[(i+1)*nx+(j  )].x-pcr[(i  )*nx+(j  )].x)/2.+.5);
            // shift of y coordinate for one vertical mesh step
            ydy = (int)((pcr[(i+1)*nx+(j+1)].y-pcr[(i  )*nx+(j+1)].y+
                        pcr[(i+1)*nx+(j  )].y-pcr[(i  )*nx+(j  )].y)/2.+.5);
            sb = sw = nb = nw = 0;
            // black at diagonal shifts
            if ((cc = ownIPL_ADJ_GetLocalColor(im,xs,ys,x-xdx-xdy,y-ydx-ydy))>=0)
            {
              sb += cc;
              nb++;
            }
            if ((cc = ownIPL_ADJ_GetLocalColor(im,xs,ys,x-xdx+xdy,y-ydx+ydy))>=0)
            {
              sb += cc;
              nb++;
            }
            if ((cc = ownIPL_ADJ_GetLocalColor(im,xs,ys,x+xdx-xdy,y+ydx-ydy))>=0)
            {
              sb += cc;
              nb++;
            }
            if ((cc = ownIPL_ADJ_GetLocalColor(im,xs,ys,x+xdx+xdy,y+ydx+ydy))>=0)
            {
              sb += cc;
              nb++;
            }
            // white at horizontal and vertical shifts
            if ((cc = ownIPL_ADJ_GetLocalColor(im,xs,ys,x-xdx,y-ydx))>=0)
            {
              sw += cc;
              nw++;
            }
            if ((cc = ownIPL_ADJ_GetLocalColor(im,xs,ys,x+xdx,y+ydx))>=0)
            {
              sw += cc;
              nw++;
            }
            if ((cc = ownIPL_ADJ_GetLocalColor(im,xs,ys,x-xdy,y-ydy))>=0)
            {
              sw += cc;
              nw++;
            }
            if ((cc = ownIPL_ADJ_GetLocalColor(im,xs,ys,x+xdy,y+ydy))>=0)
            {
              sw += cc;
              nw++;
            }
            // we should have at least two references for each color
            if ((nb>=2)&&(nw>=2))
              // get our color
              if ((cc = ownIPL_ADJ_GetLocalColor(im,xs,ys,x,y))>=0)
                // we should fit to zone closer than 1/3 from white to black
                if (3*cc>(int)((2*sw/nw)+(sb/nb)))
                { // Yes! We have found the marker
                  // save it's position, free resources, and exit
                  Marker->x = x;
                  Marker->y = y;
                  free(pcr);
                  return ERR_OK;
                }
            }
          }
    // all points searched. no marker found
    free(pcr);
    return ERR_IPL_NO_OBJECTS_FOUND;
}

// determine the central white blob automatically
MIA_RESULT_CODE IPL_ADJ_FindCentralBlob(
  SDot* Marker,               // OUT: central dot position
  const unsigned char* im,    // IN:  source image
  int xs,           // IN:  image size
  int ys)           //
{
  MIA_RESULT_CODE res = ERR_GEN_NOMEMORY;
  unsigned char *bufs[4];
  int crossnum,orderedCN,nx,ny;
  SDot *pCrosses = NULL,*orderedCrosses = NULL;

    if ((bufs[0] = (unsigned char*)malloc(xs*ys*4))==NULL)
      goto hahaxit;
    bufs[1] = bufs[0]+xs*ys;
    bufs[2] = bufs[1]+xs*ys;
    bufs[3] = bufs[2]+xs*ys;
    if ((pCrosses = (SDot*)malloc(sizeof(SDot)*(crossnum = 200)))==NULL)
      goto hahaxit;
    res = IPL_ADJ_GetChessCrosses2(
      pCrosses,           // OUT: crosses found
      &crossnum,       // IN/OUT: reserved space/number of actually found crosses
      im,  // IN: source image
      xs,         // IN: image width in pixels
      ys,         // IN: image height in pixels
      bufs);    // temporary buffers
    if (res)
      goto hahaxit;
    res = IPL_ADJ_RefineChessCrosses(
      pCrosses,           // OUT: crosses
      crossnum,        // IN: number of crosses
      im,  // IN: source image
      xs,         // IN: image width in pixels
      ys,
      0, 0);         // IN: window size (0 means function guesses itself)
    if (res)
      goto hahaxit;
    orderedCN = crossnum;
    res = IPL_ADJ_OrderPointsInMesh(
      &orderedCrosses,    // OUT: outputted points in a raster (rectangle mesh) form
      pCrosses,  // IN:  input points array
      &orderedCN,   // IN:  number of points, OUT: number of approved points
      &nx,    // OUT: horizontal size of mesh
      &ny,    // OUT: vertical size of mech
      xs,     // IN:  image width
      ys,    // IN:  image height
      NULL);
    if (res)
      goto hahaxit;
    res = IPL_ADJ_FindWhiteCircleMarker(
      Marker,             // OUT: marker position
      orderedCrosses,      // IN:  crosses array
      nx,         // IN:  crosses array size
      ny,         //
      im,  // IN:  image
      xs,         // IN:  image size
      ys);        //
hahaxit:
    if (bufs[0])
      free(bufs[0]);
    if (pCrosses)
      free(pCrosses);
    if (orderedCrosses)
      free(orderedCrosses);
    return res;
}

// calculate the value showing deviation of mesh from ideally parallelogramic
MIA_RESULT_CODE IPL_ADJ_EvaluateParallelogramicMeshQuality(
  double* quality,        // OUT: evaluated quality
  double* uncertainty,    // OUT: uncertainty (due to holes in mesh)
  const SDot* mesh,       // IN:  mesh to be evaluated
  int nx,       // IN:  dimensions of the mesh
  int ny)       // IN:
{
  SMeshLineParameters* smlp;
  double dummy_uncertainty,sum,sum_u;
  int i,*pidx;

    // check arguments
    if ((quality==NULL)||(mesh==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((nx==0)||(ny==0))
      return ERR_GEN_NO_DATA;
    if (uncertainty==NULL)
      uncertainty = &dummy_uncertainty;
    // allocate
    if ((smlp = (SMeshLineParameters*)malloc(sizeof(SMeshLineParameters)*(nx+ny)))==NULL)
      return ERR_GEN_NOMEMORY;
    if ((pidx = (int*)malloc(sizeof(int)*(nx+ny)))==NULL)
    {
      free(smlp);
      return ERR_GEN_NOMEMORY;
    }
    // vertical lines
    for (i=0;i<nx;i++)
      IPL_ADJ_EvaluateLineQuality(
        smlp+i, // OUT: calculated parameters
        mesh+i,         // IN:  start of line in mesh
        nx,        // IN:  stride between line elements
        ny,          // IN:  length of line
        pidx);      // TMP: index buffer of n elements
    // horizontal lines
    for (i=0;i<ny;i++)
      IPL_ADJ_EvaluateLineQuality(
        smlp+nx+i, // OUT: calculated parameters
        mesh+i*nx,         // IN:  start of line in mesh
        1,        // IN:  stride between line elements
        nx,          // IN:  length of line
        pidx);      // TMP: index buffer of n elements
    free(pidx);
    for (sum_u = sum = 0.,i=0;i<nx+ny;i++)
    {
      sum += smlp[i].quality;
      sum_u += smlp[i].uncertainty;
    }
    free(smlp);
    *quality = sum/(nx+ny);
    *uncertainty = sum_u/(nx+ny);
    return ERR_OK;
}

// calculate mean and stddev of chessboard cell pixel size
MIA_RESULT_CODE IPL_ADJ_CalculatePixDist(
  double* pdist_x,    // OUT: pixel size along x
  double* pdist_y,    // OUT: pixel size along y
  double* pdev,       // OUT: deviation
  const SDot* dots,   // IN:  array of dots (unsorted, unaligned)
  int num,  // IN:  size of array
  int xs,   // IN:  image size
  int ys)   // IN:
{
  SDot* ordered;
  int nx,ny,i,j,num_x,num_y;
  MIA_RESULT_CODE res;
  double d,dist_x,dist_y,dev_x,dev_y;

    // check argumats
    if (dots==NULL)
      return ERR_GEN_NULLPOINTER;
    // order points to square grid
    if ((res = IPL_ADJ_OrderPointsInMesh(
      &ordered,    // OUT: outputted points in a raster (rectangle mesh) form
      dots,  // IN:  input points array
      &num,   // IN:  number of points, OUT: number of approved points
      &nx,    // OUT: horizontal size of mesh
      &ny,    // OUT: vertical size of mech
      xs,     // IN:  image width
      ys, NULL))!=ERR_OK)
      return res;
    // calculate requred values
    num_x = num_y = 0;
    dist_x = dist_y = dev_x = dev_y = 0.;
    for (i=0;i<ny-1;i++)
      for (j=0;j<nx-1;j++)
      {
        // vertical neighbor
        if ((ordered[(i  )*nx+j].isDefined)&&
            (ordered[(i+1)*nx+j].isDefined))
        {
          d = sqrt((ordered[(i+1)*nx+j].x-ordered[(i  )*nx+j].x)*
                   (ordered[(i+1)*nx+j].x-ordered[(i  )*nx+j].x)+
                   (ordered[(i+1)*nx+j].y-ordered[(i  )*nx+j].y)*
                   (ordered[(i+1)*nx+j].y-ordered[(i  )*nx+j].y));
          num_y++;
          dist_y += d;
          dev_y += d*d;
        }
        // horizontal neighbor
        if ((ordered[(i  )*nx+(j  )].isDefined)&&
            (ordered[(i  )*nx+(j+1)].isDefined))
        {
          d = sqrt((ordered[(i  )*nx+(j+1)].x-ordered[(i  )*nx+j].x)*
                   (ordered[(i  )*nx+(j+1)].x-ordered[(i  )*nx+j].x)+
                   (ordered[(i  )*nx+(j+1)].y-ordered[(i  )*nx+j].y)*
                   (ordered[(i  )*nx+(j+1)].y-ordered[(i  )*nx+j].y));
          num_x++;
          dist_x += d;
          dev_x += d*d;
        }
      }
    // delete array
    free(ordered);
    // are there any good neighbours?
    if ((!num_x)||(!num_y))
      return ERR_GEN_NOT_INITIALISED;
    // output result
    if (pdist_x)
      *pdist_x = dist_x/num_x;
    if (pdist_y)
      *pdist_y = dist_y/num_y;
    if (pdev)
      *pdev = ((dev_x/num_x-(dist_x/num_x)*(dist_x/num_x))+
               (dev_x/num_x-(dist_x/num_x)*(dist_x/num_x)) )/2.;
    return ERR_OK;
}

// combine transformations given by XY-maps
// C[x] = (B * A)[x] = B[A[x]]   (A is done first)
MIA_RESULT_CODE IPL_ADJ_CombineTransfomations(
  double* C_x,          // OUT: destination X map
  double* C_y,          // OUT: destination Y map
  const double* A_x,    // IN:  first source X
  const double* A_y,    // IN:  first source Y
  const double* B_x,    // IN:  second source X
  const double* B_y,    // IN:  second source Y
  int xs,     // IN:  transformation map size
  int ys)     //
{
  double p,q,S,u,v,SX,SY,pq;
  int i,j,x,y;

    // check arguments
    if ((xs<2)||(ys<2))
      return ERR_GEN_NO_DATA;
    if ((C_x==NULL)||(C_y==NULL)||
        (A_x==NULL)||(A_y==NULL)||
        (B_x==NULL)||(B_y==NULL) )
      return ERR_GEN_NULLPOINTER;
    //
    for (i=0;i<(int)(ys);i++)
      for (j=0;j<(int)(xs);j++)
      {
        // calculate partitions
        if ((S = B_x[i*xs+j])>=0)
          p = S-(x = (int)(S));
        else
           p = x = -1; // just sign bad // p = S-(x = -1-int(-S));
        if ((S = B_y[i*xs+j])>=0)
          q = S-(y = (int)(S));
        else
          q = y = -1; // just sign bad // q = S-(y = -1-int(-S));
        // is B transformation item in borders?
        if ( ((x>=0)&&((x<(int)(xs))||((x==(int)(xs))&&(p==0)))) &&
             ((y>=0)&&((y<(int)(ys))||((y==(int)(ys))&&(q==0)))) )
        { // B step pixel is OK. enumerate BLI neighbors
          SX = SY = S = 0.;
          // 0,0
          if ((pq = (1.-p)*(1.-q))>0)
          {
            u = A_x[(y  )*xs+(x  )];
            v = A_y[(y  )*xs+(x  )];
            if ((u>=0)&&(u<=xs)&&(v>=0)&&(v<=ys))
            {
              SX += u*pq;
              SY += v*pq;
              S  += pq;
            }
          }
          // 0,1
          if (((pq = p*(1.-q))>0) && ((x+1)<(int)xs))
          {
            u = A_x[(y  )*xs+(x+1)];
            v = A_y[(y  )*xs+(x+1)];
            if ((u>=0)&&(u<=xs)&&(v>=0)&&(v<=ys))
            {
              SX += u*pq;
              SY += v*pq;
              S  += pq;
            }
          }
          // 1,0
          if (((pq = (1.-p)*q)>0) && ((y+1)<(int)ys))
          {
            u = A_x[(y+1)*xs+(x  )];
            v = A_y[(y+1)*xs+(x  )];
            if ((u>=0)&&(u<=xs)&&(v>=0)&&(v<=ys))
            {
              SX += u*pq;
              SY += v*pq;
              S  += pq;
            }
          }
          // 1,1
          if (((pq = p*q)>0) && ((x+1)<(int)xs) && ((y+1)<(int)ys))
          {
            u = A_x[(y+1)*xs+(x+1)];
            v = A_y[(y+1)*xs+(x+1)];
            if ((u>=0)&&(u<=xs)&&(v>=0)&&(v<=ys))
            {
              SX += u*pq;
              SY += v*pq;
              S  += pq;
            }
          }
          // save now
          if (S)
          {
            C_x[i*xs+j] = SX/S;
            C_y[i*xs+j] = SY/S;
          }
          else
          { // all A's are outbound. just assign Ci=Ai
//            C_x[i*xs+j] = A_x[y*xs+x];
//            C_y[i*xs+j] = A_y[y*xs+x];
            C_x[i*xs+j] = -1000.;
            C_y[i*xs+j] = -1000.;
          }
        }
        else
        { // outbound from array on step B. just assign Ci=Bi
//          C_x[i*xs+j] = B_x[i*xs+j];
//          C_y[i*xs+j] = B_y[i*xs+j];
          C_x[i*xs+j] = -1000.;
          C_y[i*xs+j] = -1000.;
        }
      }
    return ERR_OK;
}
