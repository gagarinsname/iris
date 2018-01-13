/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  16 March 2006                                          */
/*      Revision: 1.0.00                                                 */
/*      Purpose:  pumps - i.e. image transformation, where destination   */
/*                point value is made out of single source point value   */
/*                or each source point influences only one destination   */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/
#define __FILENUM__ 113 // __FILENUM__TAG113

#include <malloc.h>
#include <string.h>
#include "ipl.h"

// make transformation of 8x8 square to one pixel
// even interlace half-frame is taken, so 8*4 points are actually used
void IPL_PMP_PyramidMidOf8x4IL(
  unsigned char* dst,
  const unsigned char* src,
  int xs,   // source size
  int ys)
{
  int i,j;

    for (i=ys/8-1;i>=0;i--)
      for (j=xs/8-1;j>=0;j--)
        dst[i*(xs/8)+j] = (unsigned char)((
          src[(i*8  )*xs+(j*8  )]+src[(i*8+2)*xs+(j*8  )]+src[(i*8+4)*xs+(j*8  )]+src[(i*8+6)*xs+(j*8  )]+
          src[(i*8  )*xs+(j*8+1)]+src[(i*8+2)*xs+(j*8+1)]+src[(i*8+4)*xs+(j*8+1)]+src[(i*8+6)*xs+(j*8+1)]+
          src[(i*8  )*xs+(j*8+2)]+src[(i*8+2)*xs+(j*8+2)]+src[(i*8+4)*xs+(j*8+2)]+src[(i*8+6)*xs+(j*8+2)]+
          src[(i*8  )*xs+(j*8+3)]+src[(i*8+2)*xs+(j*8+3)]+src[(i*8+4)*xs+(j*8+3)]+src[(i*8+6)*xs+(j*8+3)]+
          src[(i*8  )*xs+(j*8+4)]+src[(i*8+2)*xs+(j*8+4)]+src[(i*8+4)*xs+(j*8+4)]+src[(i*8+6)*xs+(j*8+4)]+
          src[(i*8  )*xs+(j*8+5)]+src[(i*8+2)*xs+(j*8+5)]+src[(i*8+4)*xs+(j*8+5)]+src[(i*8+6)*xs+(j*8+5)]+
          src[(i*8  )*xs+(j*8+6)]+src[(i*8+2)*xs+(j*8+6)]+src[(i*8+4)*xs+(j*8+6)]+src[(i*8+6)*xs+(j*8+6)]+
          src[(i*8  )*xs+(j*8+7)]+src[(i*8+2)*xs+(j*8+7)]+src[(i*8+4)*xs+(j*8+7)]+src[(i*8+6)*xs+(j*8+7)]+ 15)/32);
}

// make transformation of 4x4 square to one pixel
// even interlace half-frame is taken, so 4*2 points are actually used
void IPL_PMP_PyramidMidOf4x2IL(
  unsigned char* dst,
  const unsigned char* src,
  int xs,   // source size
  int ys)
{
  int i,j;

    for (i=ys/4-1;i>=0;i--)
      for (j=xs/4-1;j>=0;j--)
        dst[i*(xs/4)+j] = (unsigned char)((
          src[(i*4  )*xs+(j*4  )]+src[(i*4+2)*xs+(j*4  )]+
          src[(i*4  )*xs+(j*4+1)]+src[(i*4+2)*xs+(j*4+1)]+
          src[(i*4  )*xs+(j*4+2)]+src[(i*4+2)*xs+(j*4+2)]+
          src[(i*4  )*xs+(j*4+3)]+src[(i*4+2)*xs+(j*4+3)]+ 3)/8);
}

// make transformation of 4x4 square to one pixel
// even interlace half-frame is taken, so 4*2 points are actually used
void IPL_PMP_PyramidMidOf4x2IL_int(
  int* dst,
  const unsigned char* src,
  int xs,   // source size
  int ys)
{
  int i,j;

    for (i=ys/4-1;i>=0;i--)
      for (j=xs/4-1;j>=0;j--)
        dst[i*(xs/4)+j] = (int)
          src[(i*4  )*xs+(j*4  )]+src[(i*4+2)*xs+(j*4  )]+
          src[(i*4  )*xs+(j*4+1)]+src[(i*4+2)*xs+(j*4+1)]+
          src[(i*4  )*xs+(j*4+2)]+src[(i*4+2)*xs+(j*4+2)]+
          src[(i*4  )*xs+(j*4+3)]+src[(i*4+2)*xs+(j*4+3)];
}

// make transformation of 4x4 square to one pixel
// even interlace half-frame is taken, so 4*2 points are actually used
void IPL_PMP_PyramidMidOf4x2IL_Ext(
  unsigned char* dst,       // OUT: destination
  int dststr,               // IN:  dest stride in bytes
  int x,                    // IN:  window in source
  int y,                    //
  int w,                    //
  int h,                    //
  const unsigned char* src, // IN:  source image
  int srcstr,               // IN:  source stride in bytes
  int xs,                   // IN:  source image size
  int ys)                   //
{
  int i,j;

    w = (w+3)&(~3);
    h = (h+3)&(~3);
    for (i=y;i<y+4*h;i+=4)
      for (j=x;j<x+4*w;j+=4)
        if ((i<0)||(i+3>ys)||(j<0)||(j+3>xs))
          dst[((i-y)/4)*dststr+((j-x)/4)] = 0;
        else
          dst[((i-y)/4)*dststr+((j-x)/4)] = (unsigned char)((
            src[(i  )*srcstr+(j  )]+src[(i+2)*srcstr+(j  )]+
            src[(i  )*srcstr+(j+1)]+src[(i+2)*srcstr+(j+1)]+
            src[(i  )*srcstr+(j+2)]+src[(i+2)*srcstr+(j+2)]+
            src[(i  )*srcstr+(j+3)]+src[(i+2)*srcstr+(j+3)]+ 3)/8);
}

// make transformation of 2x2 square to one pixel
// even interlace half-frame is taken, so 2*1 points are actually used
void IPL_PMP_PyramidMidOf2x1IL(
  unsigned char* dst,
  const unsigned char* src,
  int xs,   // source size
  int ys)
{
  int i,j;

    for (i=ys/2-1;i>=0;i--)
      for (j=xs/2-1;j>=0;j--)
        dst[i*(xs/2)+j] = (unsigned char)((
          src[(i*2  )*xs+(j*2  )]+
          src[(i*2  )*xs+(j*2+1)] )/2);
}

void IPL_PMP_PyramidMidOf2x2(
  unsigned char* dst,
  const unsigned char* src,
  int xs,   // source size
  int ys)
{
  int i,j;

    for (i=ys/2-1;i>=0;i--)
      for (j=xs/2-1;j>=0;j--)
        dst[i*(xs/2)+j] = (unsigned char)((
          src[(i*2  )*xs+(j*2  )]+src[(i*2+1)*xs+(j*2  )]+
          src[(i*2  )*xs+(j*2+1)]+src[(i*2+1)*xs+(j*2+1)] + 2)/4);
}

// copy odd rows to even rows / or vice versa
// BpP!=1 can be used
MIA_RESULT_CODE IPL_PMP_CopyImageStrings(
  unsigned char* dst,           // destination
  const unsigned char* src,     // source
  int xs,                       // image size
  int ys,                       //
  int BpP,                      // bytes per pixel
  int bInterpolateOdd)          // if !=0, calculate odd from even
{
  int i;

    if ((ys==1)&&(!bInterpolateOdd))
      bInterpolateOdd = 1;
    if (!bInterpolateOdd)
    { // even rows are replaced
      // zero-th row copied from first
      MIA_memcpy(dst,src+xs*BpP,xs*BpP);
      // adjust pointers as is it were even rows
      dst += xs*BpP;
      src += xs*BpP;
      ys--;
      bInterpolateOdd = 1;
    }
    if (!(ys&1))
      // last row is odd
      MIA_memcpy(dst+(ys-1)*xs*BpP,src+(ys-2)*xs*BpP,xs*BpP);
    // body.
    // copy even to even
    for (i=0;i<ys;i+=2)
      MIA_memcpy(dst+i*xs*BpP,src+i*xs*BpP,xs*BpP);
    // copy even to odd
    for (i=0;i<ys-2;i+=2)
      MIA_memcpy(dst+(i+1)*xs*BpP,src+i*xs*BpP,xs*BpP);
    return ERR_OK;
}

void IPL_PMP_MakeIntegralImage(
  int* intim,               // OUT: integral image
  const unsigned char* im,  // IN:  source image
  int xs,                   // IN:  size
  int ys)
{
  int i,j;

    intim[0] = im[0];
    for (j=1;j<xs;j++)
      intim[j] = im[j]+intim[j-1];
    for (i=1;i<ys;i++)
      intim[i*xs] = im[i*xs]+intim[(i-1)*xs];
    for (i=1;i<ys;i++)
      for (j=1;j<xs;j++)
        intim[i*xs+j] =  im   [i*xs+ j   ]    -intim[(i-1)*xs+(j-1)]
                        +intim[(i-1)*xs+ j   ]+intim[ i   *xs+(j-1)];
}

void IPL_PMP_MakeIntegralImage_prim(
  int* intim,               // OUT: integral image, stride is (xs+1) elements
  const unsigned char* im,  // IN:  source image
  int xs,                   // IN:  size
  int ys)
{
  int i,j;

    memset(intim,0,sizeof(*intim)*(xs+1));
    for (i=1;i<=ys;i++)
      intim[i*(xs+1)] = 0;
    for (i=1;i<=ys;i++)
      for (j=1;j<=xs;j++)
        intim[i*(xs+1)+j] =  im   [(i-1)* xs   +(j-1)]-intim[(i-1)*(xs+1)+(j-1)]
                            +intim[(i-1)*(xs+1)+ j   ]+intim[ i   *(xs+1)+(j-1)];
}

// summing downscaling by 2
MIA_RESULT_CODE IPL_PMP_Squeeze2_Good(
  unsigned char* dst,
  int dststr,
  const unsigned char* src,
  int srcstr,
  int xs,
  int ys,
  int BpP)
{
  int i,j,k;

    for (i=0;i<ys/2;i++)
    {
      for (j=0;j<xs/2;j++)
        for (k=0;k<BpP;k++)
          dst[j*BpP+k] = (unsigned char)(
            ((unsigned int)src[(j*2         )*BpP+k]+
                           src[(j*2+1       )*BpP+k]+
                           src[(j*2+srcstr  )*BpP+k]+
                           src[(j*2+srcstr+1)*BpP+k]+1)/4);
      dst += dststr;
      src += srcstr*2;
    }
    return ERR_OK;
}

// simple downscaling by 2
MIA_RESULT_CODE IPL_PMP_Squeeze2_Simple(
  unsigned char* dst,
  int dststr,
  const unsigned char* src,
  int srcstr,
  int xs,
  int ys,
  int BpP)
{
  int i,j,k;

    for (i=0;i<ys/2;i++)
    {
      for (j=0;j<xs/2;j++)
        for (k=0;k<BpP;k++)
          dst[j*BpP+k] = src[(j*2         )*BpP+k];
      dst += dststr;
      src += srcstr*2;
    }
    return ERR_OK;
}

// simple upscaling by 2
MIA_RESULT_CODE IPL_PMP_Expand2_Simple(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  int BpP)
{
  int i,j,k;

    for (i=0;i<ys;i++)
      for (j=0;j<xs;j++)
        for (k=0;k<BpP;k++)
          dst[((i*2  )*(xs*2)+(j*2  ))*BpP+k] =
          dst[((i*2  )*(xs*2)+(j*2+1))*BpP+k] =
          dst[((i*2+1)*(xs*2)+(j*2  ))*BpP+k] =
          dst[((i*2+1)*(xs*2)+(j*2+1))*BpP+k] =
            src[((i)*(xs)+(j))*BpP+k];
    return ERR_OK;
}

// interpolating upscaling by 2
MIA_RESULT_CODE IPL_PMP_Expand2_Good(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  int BpP)
{
  int i,j,k;

    // body
    for (i=0;i<ys-1;i++)
      for (j=0;j<xs-1;j++)
        for (k=0;k<BpP;k++)
        {
          dst[((i*2  )*(xs*2)+(j*2  ))*BpP+k] = src[((i  )*(xs)+(j  ))*BpP+k];
          dst[((i*2  )*(xs*2)+(j*2+1))*BpP+k] =
            (unsigned char)(
            ((unsigned int)src[((i  )*(xs)+(j  ))*BpP+k]+
                           src[((i  )*(xs)+(j+1))*BpP+k] )/2);
          dst[((i*2+1)*(xs*2)+(j*2  ))*BpP+k] =
            (unsigned char)(
            ((unsigned int)src[((i  )*(xs)+(j  ))*BpP+k]+
                           src[((i+1)*(xs)+(j  ))*BpP+k] )/2);
          dst[((i*2+1)*(xs*2)+(j*2+1))*BpP+k] =
            (unsigned char)(
            ((unsigned int)src[((i  )*(xs)+(j  ))*BpP+k]+
                           src[((i  )*(xs)+(j+1))*BpP+k]+
                           src[((i+1)*(xs)+(j  ))*BpP+k]+
                           src[((i+1)*(xs)+(j+1))*BpP+k]+1)/4);
        }
    // bottom
    i = ys-1;
    for (j=0;j<xs-1;j++)
      for (k=0;k<BpP;k++)
      {
        dst[((i*2  )*(xs*2)+(j*2  ))*BpP+k] =
        dst[((i*2+1)*(xs*2)+(j*2  ))*BpP+k] =
          src[((i  )*(xs)+(j  ))*BpP+k];
        dst[((i*2  )*(xs*2)+(j*2+1))*BpP+k] =
        dst[((i*2+1)*(xs*2)+(j*2+1))*BpP+k] =
          (unsigned char)(
          ((unsigned int)src[((i  )*(xs)+(j  ))*BpP+k]+
                         src[((i  )*(xs)+(j+1))*BpP+k] )/2);
      }
    // right border
    j = xs-1;
    for (i=0;i<ys-1;i++)
      for (k=0;k<BpP;k++)
      {
        dst[((i*2  )*(xs*2)+(j*2  ))*BpP+k] =
        dst[((i*2  )*(xs*2)+(j*2+1))*BpP+k] =
          src[((i  )*(xs)+(j  ))*BpP+k];
        dst[((i*2+1)*(xs*2)+(j*2  ))*BpP+k] =
        dst[((i*2+1)*(xs*2)+(j*2+1))*BpP+k] =
          (unsigned char)(
          ((unsigned int)src[((i  )*(xs)+(j  ))*BpP+k]+
                         src[((i+1)*(xs)+(j  ))*BpP+k] )/2);
      }
    // bottom-right corner
    i = ys-1;
    for (k=0;k<BpP;k++)
      dst[((i*2  )*(xs*2)+(j*2  ))*BpP+k] =
      dst[((i*2  )*(xs*2)+(j*2+1))*BpP+k] =
      dst[((i*2+1)*(xs*2)+(j*2  ))*BpP+k] =
      dst[((i*2+1)*(xs*2)+(j*2+1))*BpP+k] =
        src[((i  )*(xs)+(j  ))*BpP+k];
    return ERR_OK;
}

// combining half-frames
MIA_RESULT_CODE IPL_PMP_CombineHalfFrames(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  int bOddIsFirst)
{
  int i;

    // check arguments
    if ((dst==NULL)||(src==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((xs<=0)||(ys<=0))
      return ERR_GEN_NO_DATA;
    // decide the copy order
    if (bOddIsFirst)
    { // copy first half to odd strings, second half to even strings
      for (i=0;i<(ys+1)/2;i++)
        MIA_memcpy(dst+(2*i+1)*xs,src+ i         *xs,xs);
      for (i=0;i<(ys  )/2;i++)
        MIA_memcpy(dst+ 2*i   *xs,src+(i+(ys+1)/2)*xs,xs);
    }
    else
    { // copy first half to even strings, second half to odd strings
      for (i=0;i<(ys+1)/2;i++)
        MIA_memcpy(dst+ 2*i   *xs,src+ i         *xs,xs);
      for (i=0;i<(ys  )/2;i++)
        MIA_memcpy(dst+(2*i+1)*xs,src+(i+(ys+1)/2)*xs,xs);
    }
    return ERR_OK;
}

// splitting half-frames
// order 0-1-2-3-4-5 changes to 0-2-4-1-3-5
MIA_RESULT_CODE IPL_PMP_SplitHalfFrames(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  int bOddIsFirst)
{
  int i;

    // check arguments
    if ((dst==NULL)||(src==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((xs<=0)||(ys<=0))
      return ERR_GEN_NO_DATA;
    // decide the copy order
    if (bOddIsFirst)
    { // copy odd strings to the first half, even to the second
      for (i=0;i<(ys+1)/2;i++)
        MIA_memcpy(dst+ i          *xs,src+(2*i+1)*xs,xs);
      for (i=0;i<(ys  )/2;i++)
        MIA_memcpy(dst+(i+(ys+1)/2)*xs,src+ 2*i   *xs,xs);
    }
    else
    { // copy even strings to the first half, odd to the second
      for (i=0;i<(ys+1)/2;i++)
        MIA_memcpy(dst+ i          *xs,src+ 2*i   *xs,xs);
      for (i=0;i<(ys  )/2;i++)
        MIA_memcpy(dst+(i+(ys+1)/2)*xs,src+(2*i+1)*xs,xs);
    }
    return ERR_OK;
}

// interpolate odd rows with even rows / or vice versa
// BpP!=1 can be used
MIA_RESULT_CODE IPL_PMP_Interpolate(
  unsigned char* dst,           // destination
  const unsigned char* src,     // source
  int xs,                       // image size
  int ys,                       //
  int BpP,                      // bytes per pixel
  int bInterpolateOdd)          // odd = f(even)
{
  int i,j,k;

    if (src==NULL)
      src = dst;
    if ((ys==1)&&(!bInterpolateOdd))
      bInterpolateOdd = 1;
    if (!bInterpolateOdd)
    { // even rows are replaced
      // zero-th row copied from first
      MIA_memcpy(dst,src+xs*BpP,xs*BpP);
      // adjust pointers as is it were even rows
      dst += xs*BpP;
      src += xs*BpP;
      ys--;
      bInterpolateOdd = 1;
    }
    if (!(ys&1))
      // last row is odd
      MIA_memcpy(dst+(ys-1)*xs*BpP,src+(ys-2)*xs*BpP,xs*BpP);
    // copy even to even
    for (i=0;i<ys;i+=2)
      MIA_memcpy(dst+i*xs*BpP,src+i*xs*BpP,xs*BpP);
    // interpolate even to odd
    for (i=0;i<ys-2;i+=2)
      for (j=0;j<xs;j++)
        for (k=0;k<BpP;k++)
          dst[((i+1)*xs+(j))*BpP+k] = (unsigned char)((
               ((unsigned int)(src[((i  )*xs+(j))*BpP+k]))+
               ((unsigned int)(src[((i+2)*xs+(j))*BpP+k])) )/2);
    return ERR_OK;
}

// crop piece of source image to destination image
void IPL_PMP_CropCentering(
  unsigned char* dst,         // IN:  destination image
  int dxs,                    // IN:  dst image size
  int dys,                    //
  const unsigned char* src,   // IN:  source image
  int sxs,                    // IN:  src image size
  int sys,                    //
  int xc,                     // IN:  center coordinates
  int yc)                     //
{
  int x,y,i;

    if ((dxs>sxs)||(dys>sys))
      return;
    x = xc-dxs/2;
    y = yc-dys/2;
    if (x<0)
      x = 0;
    if (x+dxs>sxs)
      x = sxs-dxs;
    if (y<0)
      y = 0;
    if (y+dys>sys)
      y = sys-dys;
    for (i=0;i<dys;i++)
      memcpy(dst+i*dxs,src+(y+i)*sxs+x,dxs);
}

// crop piece of source image to destination image
void IPL_PMP_Crop(
  unsigned char* dst,         // IN:  destination image
  int dstr,                   // IN:  dst image stride
  int dxs,                    // IN:  dst image size
  int dys,                    //
  int xplace,                 // IN:  position to insert the ROI
  int yplace,                 //
  const unsigned char* src,   // IN:  source image
  int sstr,                   // IN:  src image stride
  int sxs,                    // IN:  src image size
  int sys,                    //
  int x,                      // IN:  source ROI
  int y,                      //
  int w,                      //
  int h)                      //
{
  int idx;

    // natural limitations
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
    // limitations of source
    if (x+w>sxs)
      w = sxs-x;
    if (y+h>sys)
      h = sys-y;
    // limitations of destination
    if (xplace+w>dxs)
      w = dxs-xplace;
    if (yplace+h>dys)
      h = dys-yplace;
    // check for adequacy
    if ((w<=0)||(h<=0))
      return;
    // crop
    for (idx=0;idx<h;idx++)
      memcpy(dst+(yplace+idx)*dstr+xplace,src+(y+idx)*sstr+x,w);
}

staticFunc void ownIPL_makemix(
  unsigned char* dst,
  double frac,
  const unsigned char* src1,
  const unsigned char* src2,
  int sz)
{
  int i;

    if (frac==0.)
      memcpy(dst,src1,sz);
    else
      for (i=0;i<sz;i++)
        dst[i] = (unsigned char)(src1[i]*(1.-frac)+src2[i]*frac+.5);
}

staticFunc void ownIPL_darken2(
  unsigned char* dst,
  char shift,
  int sz)
{
  int i,v;

    for (i=0;i<sz;i++)
    {
      v = (int)shift*2+((unsigned int)dst[i]);
      if (v<0)
        v = 0;
      if (v>255)
        v = 255;
      dst[i] = (unsigned char)v;
    }
}

// make time resampling transformation
MIA_RESULT_CODE IPL_PMP_ResampleImageSequence(
  unsigned char* pDst,        // OUT: resampled sequence
  int nFrDst,                 // IN:  number of full frames in the resampled sequence
  const char* pBrtShifts,     // IN:  brightness shifts, may be NULL
  const unsigned char* pSrc,  // IN:  source image sequence
  int nFrSrc,                 // IN:  number of full frames in source
  int nInserted,              // IN:  number of duplicated half-frames in the beg.
  int xs,                     // IN:  image width in pixels
  int ys,                     // IN:  image height in pixels
  int bOddIsFirst)            // IN:  order of half-frames
{
  MIA_RESULT_CODE res=ERR_OK;
  unsigned char* pTmp=NULL;
  int nFr,intg,nHFrSrc,nHFrDst;
  double x,frac;

    nHFrDst = 2*nFrDst;
    nHFrSrc = 2*nFrSrc;
    // check arguments
    if ((pDst==NULL)||(pSrc==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((nFrDst<=0)||(nFrSrc<=0)||(xs<=0)||(ys<=0))
      return ERR_GEN_NO_DATA;
    if ((nFrDst>5000)||(nInserted>=nHFrDst)||(ys&1))
      return ERR_GEN_INVALID_PARAMETER;
    if (((int64)xs*ys*nHFrSrc>=0x100000000LL)||((int64)xs*ys*nHFrDst>=0x100000000LL))
      return ERR_GEN_INSUFFICIENT_BUFFER;
    // allocate temp memory
    intg = (nFrSrc>nFrDst)?nFrSrc:nFrDst;
    if ((pTmp = (unsigned char*)malloc(xs*ys*intg))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto IPL_ResampleImageSequence_exit;
    }
    // split interlaced frames
    for (nFr=0;nFr<nFrSrc;nFr++)
      IPL_PMP_SplitHalfFrames(pTmp+nFr*xs*ys,pSrc+nFr*xs*ys,xs,ys,bOddIsFirst);
    // time resampling
    for (nFr=0;nFr<nHFrDst;nFr++)
    {
      // FrSrc moments -> HFrDst moments, endings should meet
      x = ((double)(nFr*(nFrSrc-1)))/(nHFrDst-1);
      intg = (int)x;
      frac = x-intg;
      // mix the data
      ownIPL_makemix(pDst+nFr*xs*ys/2,frac,
        pTmp+ intg   *xs*ys+(nFr&1)*xs*ys/2,
        pTmp+(intg+1)*xs*ys+(nFr&1)*xs*ys/2, xs*ys/2);
    }
    // insert frames
    memmove(pDst+nInserted*xs*ys/2,pDst,(nHFrDst-nInserted)*xs*ys/2);
    for (nFr=0;nFr<nInserted;nFr++)
      memcpy(pDst+nFr*xs*ys/2,pDst,xs*ys/2);
    // darken the frame
    if (pBrtShifts)
      for (nFr=0;nFr<nHFrDst;nFr++)
        if (pBrtShifts[nFr])
          ownIPL_darken2(pDst+nFr*xs*ys/2,pBrtShifts[nFr],xs*ys/2);
    // reorder data to temp buffer
    memcpy(pTmp,pDst,xs*ys*nFrDst);
    // combine interlaced
    for (nFr=0;nFr<nFrDst;nFr++)
      IPL_PMP_CombineHalfFrames(pDst+nFr*xs*ys,pTmp+nFr*xs*ys,xs,ys,bOddIsFirst);
IPL_ResampleImageSequence_exit:;
    if (pTmp)
      free(pTmp);
    return res;
}

// make time resampling transformation
MIA_RESULT_CODE IPL_PMP_ResampleImageSequence_Timings(
  unsigned char* pDst,        // OUT: resampled sequence
  int nFrDst,                 // IN:  number of full frames in the resampled sequence
  const char* pBrtShifts,     // IN:  brightness shifts, may be NULL
  const unsigned char* pSrc,  // IN:  source image sequence
  const int* pInTim,          // IN:  input sequence times
  int nFrSrc,                 // IN:  number of full frames in source
  int nInserted,              // IN:  number of duplicated half-frames in the beg.
  int xs,                     // IN:  image width in pixels
  int ys,                     // IN:  image height in pixels
  int bOddIsFirst)            // IN:  order of half-frames
{
  MIA_RESULT_CODE res=ERR_OK;
  unsigned char* pTmp=NULL;
  int nFr,nHFrSrc,nHFrDst,nLastSmaller;
  double frac,dInputInterval,dCurrentTime;

    nHFrDst = 2*nFrDst;
    nHFrSrc = 2*nFrSrc;
    // check arguments
    if ((pDst==NULL)||(pSrc==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((nFrDst<=0)||(nFrSrc<=0)||(xs<=0)||(ys<=0))
      return ERR_GEN_NO_DATA;
    if ((nFrDst>5000)||(nInserted>=nHFrDst)||(-nInserted>=nHFrDst)||(ys&1))
      return ERR_GEN_INVALID_PARAMETER;
    if (((int64)xs*ys*nHFrSrc>=0x100000000LL)||((int64)xs*ys*nHFrDst>=0x100000000LL))
      return ERR_GEN_INSUFFICIENT_BUFFER;
    // check income
    for (nFr=0;nFr<nFrSrc-1;nFr++)
      if (pInTim[nFr+1]<=pInTim[nFr])
        return ERR_GEN_BAD_DATA;  // later mark has smaller time
    // allocate temp memory
    nFr = (nFrSrc>nFrDst)?nFrSrc:nFrDst;
    if ((pTmp = (unsigned char*)malloc(xs*ys*nFr))==NULL)
    {
      res = ERR_GEN_NOMEMORY;
      goto IPL_PMP_ResampleImageSequence_Timings_exit;
    }
    // split interlaced frames
    for (nFr=0;nFr<nFrSrc;nFr++)
      IPL_PMP_SplitHalfFrames(pTmp+nFr*xs*ys,pSrc+nFr*xs*ys,xs,ys,bOddIsFirst);
    // get input interval
    dInputInterval = pInTim[nFrSrc-1]-pInTim[0];
    nLastSmaller = 0;
    // time resampling
    for (nFr=0;nFr<nHFrDst;nFr++)
    {
      // calculate input time matching output tick
      dCurrentTime = nFr*dInputInterval/(nHFrDst-1);
      // find the tick having biggest value, smaller than current input time
      for (;(pInTim[nLastSmaller]<=dCurrentTime)&&(nLastSmaller<nFrSrc);nLastSmaller++);
      nLastSmaller--;
      // calculate the share
      frac = (dCurrentTime-pInTim[nLastSmaller])/
             (pInTim[nLastSmaller+1]-pInTim[nLastSmaller]);
      // mix the data
      if (nLastSmaller<(nFrSrc-1))
        ownIPL_makemix(pDst+nFr*xs*ys/2,frac,
          pTmp+ nLastSmaller   *xs*ys+(nFr&1)*xs*ys/2,
          pTmp+(nLastSmaller+1)*xs*ys+(nFr&1)*xs*ys/2, xs*ys/2);
      else 
        memcpy(pDst+nFr*xs*ys/2, pTmp+ nLastSmaller*xs*ys+(nFr&1)*xs*ys/2, xs*ys/2);
    }
    // insert frames
    if (nInserted>0)
    {
      memmove(pDst+nInserted*xs*ys/2,pDst,(nHFrDst-nInserted)*xs*ys/2);
      for (nFr=0;nFr<nInserted;nFr++)
        memcpy(pDst+nFr*xs*ys/2,pDst,xs*ys/2);
    }
    else
      if (nInserted<0)
      {
        nInserted = -nInserted;
        memmove(pDst,pDst+nInserted*xs*ys/2,(nHFrDst-nInserted)*xs*ys/2);
        for (nFr=nHFrDst-nInserted;nFr<nHFrDst;nFr++)
          memcpy(pDst+nFr*xs*ys/2,pDst+(nHFrDst-nInserted-1)*xs*ys/2,xs*ys/2);
      }
    // darken the frame
    if (pBrtShifts)
      for (nFr=0;nFr<nHFrDst;nFr++)
        if (pBrtShifts[nFr])
          ownIPL_darken2(pDst+nFr*xs*ys/2,pBrtShifts[nFr],xs*ys/2);
    // reorder data to temp buffer
    memcpy(pTmp,pDst,xs*ys*nFrDst);
    // combine interlaced
    for (nFr=0;nFr<nFrDst;nFr++)
      IPL_PMP_CombineHalfFrames(pDst+nFr*xs*ys,pTmp+nFr*xs*ys,xs,ys,bOddIsFirst);
IPL_PMP_ResampleImageSequence_Timings_exit:;
    if (pTmp)
      free(pTmp);
    return res;
}

// image normalization by histogram stretch
MIA_RESULT_CODE IPL_PMP_Normalize(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  int lo_cut,     // low quantile cutoff, in 1/1000 shares
  int hi_cut)     // high quantile cutoff, in 1/1000 shares
{
  int ix,hist[256],nDnTh,nUpTh,sum,newBr;//,histSum;

    // check arguments
    if ((dst==NULL)||(src==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((xs<=0)||(ys<=0))
      return ERR_GEN_EMPTY_ROI;
    if ((lo_cut<0)||(lo_cut>=hi_cut)||(hi_cut>1000))
      return ERR_GEN_INVALID_PARAMETER;
    // fill in histogram
    memset(hist, 0, sizeof(hist));
    for (ix=0; ix<ys*xs; ix++)
      hist[src[ix]]++;
    // find absolute thresholds, as 1%-quantiles from top and bottom of hist
    sum = xs*ys*lo_cut;
    for (nDnTh=0;(nDnTh<256)&&(sum>0);nDnTh++)
      sum -= 1000*hist[nDnTh];
    sum = xs*ys*(1000-hi_cut);
    for (nUpTh=255;(nUpTh>=0)&&(sum>0);nUpTh--)
      sum -= 1000*hist[nUpTh];
    if (nDnTh==nUpTh)
    { // thresholds match - degenerate case
      memset(dst,xs*ys,nDnTh);
      return ERR_OK;
    }
    // normalize
    for (ix=0; ix<ys*xs; ix++)
    {
      newBr = (int)(255.*(src[ix]-nDnTh)/(nUpTh-nDnTh));
      if (newBr>255)
        newBr = 255;
      if (newBr<0)
        newBr=0;
      dst[ix] = (unsigned char)newBr;
    }
    return ERR_OK;
}

MIA_RESULT_CODE IPL_PMP_ScaleYUVtoGRAY(
  unsigned char* dst, // destination image
  int dststr,         // stride in bytes
  int W,              // destination image size
  int H,              //
  const unsigned char* src,   // source image
  int srcstr,         // stride in bytes
  int xs,             // source image size
  int ys)             //
{
  double A[4];

    if (W&3)
      return ERR_GEN_INVALID_PARAMETER;
    A[0] = (double)xs/W;
    A[3] = (double)ys/H;
    A[1] = A[2] = 0.;
    return IPL_GEOM_AffineTransformBLI_YUV422toGRAY(
      dst,         // destination image
      dststr,        // stride in bytes
      src,   // source image
      srcstr,        // stride in bytes
      xs,            // source image size
      ys,            //
      0.,                   // source ROI
      0.,                   //
      &(A[0]),            // affine transformation matrix
      W,             // destination ROI size
      H);            //
}

MIA_RESULT_CODE IPL_PMP_ScaleBAYERtoGRAY(
  unsigned char* dst, // destination image
  int dststr,         // stride in bytes
  int W,              // destination image size
  int H,              //
  const unsigned char* src,   // source image
  int srcstr,         // stride in bytes
  int xs,             // source image size
  int ys)             //
{
  double A[4];

    if (xs&1 || ys&1)
      return ERR_GEN_INVALID_PARAMETER;
    A[0] = (double)xs/W;
    A[3] = (double)ys/H;
    A[1] = A[2] = 0.;
    return IPL_GEOM_AffineTransformBLI_BAYERtoGRAY(
      dst,         // destination image
      dststr,        // stride in bytes
      src,   // source image
      srcstr,        // stride in bytes
      xs,            // source image size
      ys,            //
      0.,                   // source ROI
      0.,                   //
      &(A[0]),            // affine transformation matrix
      W,             // destination ROI size
      H);            //
}

MIA_RESULT_CODE IPL_PMP_ScaleBAYERtoGRAYcrop(
  unsigned char* dst, // destination image
  int dststr,         // stride in bytes
  int W,              // destination image size
  int H,              //
  const unsigned char* src,   // source image
  int srcstr,         // stride in bytes
  int xs,             // source image size
  int ys,
  int roix,                   // source ROI : x
  int roiy,					  //              y
  int roiw,                   // ROI width
  int roih					  //	 height
  )						      //
{
	double A[4];

    if (xs&1 || ys&1)
      return ERR_GEN_INVALID_PARAMETER;
    if (roix>=xs || roiy>=ys || roiw<=0 || roih<=0 || W<=0 || H<=0)
      return ERR_GEN_INVALID_PARAMETER;


    A[0] = (double)(xs-(xs-roiw))/W;
    A[3] = (double)(ys-(ys-roih))/H;
    A[1] = A[2] = 0.;
	//return IPL_GEOM_AffineTransformBLI_BAYERtoGRAY(
	return IPL_GEOM_AffineTransform_BAYERtoGRAYcrop(
      dst,         // destination image
      dststr,        // stride in bytes
      src+(roix+roiy*srcstr),   // source image
      srcstr,        // stride in bytes
      xs-roix,            // source image size
      ys-roiy,            //
      (double) 0.,   // source ROI
      (double) 0.,   //
      &(A[0]),       // affine transformation matrix
      W,             // destination ROI size
      H);            //
}

MIA_RESULT_CODE IPL_PMP_ScaleRGBtoGRAY(
  unsigned char* dst, // destination image
  int dststr,         // stride in bytes
  int W,              // destination image size
  int H,              //
  const unsigned char* src,   // source image
  int srcstr,         // stride in bytes
  int xs,             // source image size
  int ys)             //
{
  double A[4];

    if (W&3)
      return ERR_GEN_INVALID_PARAMETER;
    A[0] = (double)xs/W;
    A[3] = (double)ys/H;
    A[1] = A[2] = 0.;
    return IPL_GEOM_AffineTransformBLI_RGBtoGRAY(
      dst,         // destination image
      dststr,        // stride in bytes
      src,   // source image
      srcstr,        // stride in bytes
      xs,            // source image size
      ys,            //
      0.,                   // source ROI
      0.,                   //
      &(A[0]),            // affine transformation matrix
      W,             // destination ROI size
      H);            //
}

MIA_RESULT_CODE IPL_PMP_ScaleGrayToGray(
  unsigned char* dst, // destination image
  int dststr,         // stride in bytes
  int W,              // destination image size
  int H,              //
  const unsigned char* src,   // source image
  int srcstr,         // stride in bytes
  int xs,             // source image size
  int ys)             //
{
//    if (W&3)
  //    return ERR_GEN_INVALID_PARAMETER;
    return IPL_GEOM_RotateResizeBLI(
      dst,       // destination image
      dststr,               // stride in bytes
      src, // source image
      srcstr,               // stride in bytes
      xs,                   // source image size
      ys,                   //
      0.f,                  // source ROI
      0.f,                  //
      (float)xs,                  //
      (float)ys,                  //
      0.f,                  // angle of source ROI tilt
      W,                    // destination ROI size
      H);                    //
}

