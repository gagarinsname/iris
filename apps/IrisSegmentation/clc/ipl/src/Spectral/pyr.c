/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  10 November 2006                                       */
/*      Revision: 1.0.00                                                 */
/*      Purpose:  pyramidal representation functions                     */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/
#define __FILENUM__ 90 // __FILENUM__TAG90

#include <math.h>
#include "bpl.h"
#include "ipl.h"

MIA_RESULT_CODE IPL_PYR_SplitNNI(
  unsigned char* dstPyr,    // OUT: downscale2-pyramid
  unsigned char* srcDif,    // IN/OUT: source image/difference
  int xs,                   // IN: source image size
  int ys)                   //
{
  int i,j,v;

    // pairs of strings
    for (i=0;i<ys/2;i++)
    {
      // pairs of columns
      for (j=0;j<xs/2;j++)
      {
        v = srcDif[   0];
        *dstPyr++ = (unsigned char)v;
        srcDif[   0] = (unsigned char)(v-srcDif[   0]);
        srcDif[   1] = (unsigned char)(v-srcDif[   1]);
        srcDif[xs  ] = (unsigned char)(v-srcDif[xs  ]);
        srcDif[xs+1] = (unsigned char)(v-srcDif[xs+1]);
        srcDif += 2;
      }
      // one possible unpaired column
      if (xs&1)
      {
        v = srcDif[   0];
        *dstPyr++ = (unsigned char)v;
        srcDif[   0] = (unsigned char)(v-srcDif[   0]);
        srcDif[xs  ] = (unsigned char)(v-srcDif[xs  ]);
        srcDif += 1;
      }
      srcDif += xs;
    }
    // one possible unpaired string
    if (ys&1)
    {
      // pairs of columns
      for (j=0;j<xs/2;j++)
      {
        v = srcDif[   0];
        *dstPyr++ = (unsigned char)v;
        srcDif[   0] = (unsigned char)(v-srcDif[   0]);
        srcDif[   1] = (unsigned char)(v-srcDif[   1]);
        srcDif += 2;
      }
      // one possible unpaired column
      if (xs&1)
      {
        *dstPyr++ = srcDif[   0];
        srcDif[   0] = 0;
        srcDif += 1;
      }
    }
    return ERR_OK;
}

MIA_RESULT_CODE IPL_PYR_SplitBLI(
  unsigned char* dstPyr,    // OUT: downscale2-pyramid
  unsigned char* srcDif,    // IN/OUT: source image/difference
  int x,                   // IN: source image size
  int y)                   //
{
  size_t i,j,v;
  size_t xs = x;
  size_t ys = y;

    // pairs of strings
    for (i=0;i<ys/2;i++)
    {
      // pairs of columns
      for (j=0;j<xs/2;j++)
      {
        v = ((size_t)(srcDif[   0])+(size_t)(srcDif[   1])+
             (size_t)(srcDif[xs  ])+(size_t)(srcDif[xs+1])+2)/4;
        *dstPyr++ = (unsigned char)v;
        srcDif[   0] = (unsigned char)(v-srcDif[   0]);
        srcDif[   1] = (unsigned char)(v-srcDif[   1]);
        srcDif[xs  ] = (unsigned char)(v-srcDif[xs  ]);
        srcDif[xs+1] = (unsigned char)(v-srcDif[xs+1]);
        srcDif += 2;
      }
      // one possible unpaired column
      if (xs&1)
      {
        v = ((size_t)(srcDif[   0])+(size_t)(srcDif[xs  ])+1)/2;
        *dstPyr++ = (unsigned char)v;
        srcDif[   0] = (unsigned char)(v-srcDif[   0]);
        srcDif[xs  ] = (unsigned char)(v-srcDif[xs  ]);
        srcDif += 1;
      }
      srcDif += xs;
    }
    // one possible unpaired string
    if (ys&1)
    {
      // pairs of columns
      for (j=0;j<xs/2;j++)
      {
        v = ((size_t)(srcDif[   0])+(size_t)(srcDif[   1])+1)/2;
        *dstPyr++ = (unsigned char)v;
        srcDif[   0] = (unsigned char)(v-srcDif[   0]);
        srcDif[   1] = (unsigned char)(v-srcDif[   1]);
        srcDif += 2;
      }
      // one possible unpaired column
      if (xs&1)
      {
        *dstPyr++ = srcDif[   0];
        srcDif[   0] = 0;
        srcDif += 1;
      }
    }
    return ERR_OK;
}

int IPL_PYR_GetEnergyBLI2(
  unsigned char* dstPyr,        // OUT: downscale2-pyramid
  const unsigned char* srcDif,  // IN:  source image
  int xs,                   // IN: source image size
  int ys,                   //
  int* buf_of_512ints)
{
  int i,j,v,dx,dy,ysp,xsp;
  unsigned char* dstPyr_stor;
  const unsigned char* srcDif_stor;

    dstPyr_stor = dstPyr;
    srcDif_stor = srcDif;
    MIA_memset(buf_of_512ints,0,sizeof(buf_of_512ints[0])*512);
    // first pass - calculating pyramid
    // pairs of strings
    for (i=0;i<ys/2;i++)
    {
      // pairs of columns
      for (j=0;j<xs/2;j++)
      {
        v = ((int)(srcDif[   0])+(int)(srcDif[   1])+
             (int)(srcDif[xs  ])+(int)(srcDif[xs+1])+2)/4;
        *dstPyr++ = (unsigned char)v;
        srcDif += 2;
      }
      // one possible unpaired column
      if (xs&1)
      {
        v = ((int)(srcDif[   0])+(int)(srcDif[xs  ])+1)/2;
        *dstPyr++ = (unsigned char)v;
        srcDif += 1;
      }
      srcDif += xs;
    }
    // one possible unpaired string
    if (ys&1)
    {
      // pairs of columns
      for (j=0;j<xs/2;j++)
      {
        v = ((int)(srcDif[   0])+(int)(srcDif[   1])+1)/2;
        *dstPyr++ = (unsigned char)v;
        srcDif += 2;
      }
      // one possible unpaired column
      if (xs&1)
      {
        *dstPyr++ = srcDif[   0];
        srcDif += 1;
      }
    }
    // second pass - calculating difference histogram
    dstPyr = dstPyr_stor;
    srcDif = srcDif_stor;
    xsp = (xs+1)/2;
    ysp = (ys+1)/2;
    // pairs of strings
    for (i=0;i<ys/2;i++)
    {
      // pairs of columns
      for (j=0;j<xs/2;j++)
      {
        dx = (j       )?(   -1):0;
        dy = (i       )?(-xsp):0;
        v = (5*dstPyr[0]+dstPyr[dx+dy]+dstPyr[dx]+dstPyr[dy]+4)/8;
        buf_of_512ints[256+(v-srcDif[   0])]++;
        dx = (j<xsp-1)?(    1):0;
        dy = (i       )?(-xsp):0;
        v = (5*dstPyr[0]+dstPyr[dx+dy]+dstPyr[dx]+dstPyr[dy]+4)/8;
        buf_of_512ints[256+(v-srcDif[   1])]++;
        dx = (j       )?(   -1):0;
        dy = (i<ysp-1)?( xsp):0;
        v = (5*dstPyr[0]+dstPyr[dx+dy]+dstPyr[dx]+dstPyr[dy]+4)/8;
        buf_of_512ints[256+(v-srcDif[xs  ])]++;
        dx = (j<xsp-1)?(    1):0;
        dy = (i<ysp-1)?( xsp):0;
        v = (5*dstPyr[0]+dstPyr[dx+dy]+dstPyr[dx]+dstPyr[dy]+4)/8;
        buf_of_512ints[256+(v-srcDif[xs+1])]++;
        dstPyr++;
        srcDif += 2;
      }
      // one possible unpaired column
      if (xs&1)
      {
        dy = (i       )?(-xsp):0;
        v = (3*dstPyr[0]+dstPyr[dy]+2)/4;
        buf_of_512ints[256+(v-srcDif[   0])]++;
        dy = (i<ysp-1)?( xsp):0;
        v = (3*dstPyr[0]+dstPyr[dy]+2)/4;
        buf_of_512ints[256+(v-srcDif[xs  ])]++;
        dstPyr++;
        srcDif += 1;
      }
      srcDif += xs;
    }
    // one possible unpaired string
    if (ys&1)
    {
      // pairs of columns
      for (j=0;j<xs/2;j++)
      {
        dx = (j       )?(   -1):0;
        v = (3*dstPyr[0]+dstPyr[dx]+2)/4;
        buf_of_512ints[256+(v-srcDif[   0])]++;
        dx = (j<xsp-1)?(    1):0;
        v = (3*dstPyr[0]+dstPyr[dx]+2)/4;
        buf_of_512ints[256+(v-srcDif[   1])]++;
        dstPyr++;
        srcDif += 2;
      }
    }
    // calculate difference energy
    for (i=v=0;i<512;i++)
      v += buf_of_512ints[i]*(i-256)*(i-256);
    return (int)(sqrt(((double)v)/(xs*ys))*256);
}

// same as PYR_GetEnergyBLI2 but only inside masked area rather that in all image
int IPL_PYR_GetEnergyBLI2_withMask(
  unsigned char* dstPyr,        // OUT: downscale2-pyramid
  const unsigned char* srcDif,  // IN:  source image
  unsigned char* dstMask,       // OUT: downscale2-mask
  const unsigned char* srcMask, // IN:  source mask
  int xs,                   // IN: source image size
  int ys,                   //
  int* buf_of_512ints)
{
  int i,j,v,dx,dy,ysp,xsp,Mass=0;
  unsigned char *dstPyr_stor,*dstMask_stor;
  const unsigned char *srcDif_stor,*srcMask_stor;

    dstPyr_stor = dstPyr;
    srcDif_stor = srcDif;
    dstMask_stor = dstMask;
    srcMask_stor = srcMask;
    MIA_memset(buf_of_512ints,0,sizeof(buf_of_512ints[0])*512);
    // first pass - calculating pyramid
    // pairs of strings
    for (i=0;i<ys/2;i++)
    {
      // pairs of columns
      for (j=0;j<xs/2;j++)
      {
        // image
        v = ((int)(srcDif[   0])+(int)(srcDif[   1])+
             (int)(srcDif[xs  ])+(int)(srcDif[xs+1])+2)/4;
        *dstPyr++ = (unsigned char)v;
        srcDif += 2;
        // mask
        v = (int)(srcMask[   0])+(int)(srcMask[   1])+
            (int)(srcMask[xs  ])+(int)(srcMask[xs+1]);
//        *dstMask++ = (unsigned char)((v)?255:0);
*dstMask++ = (unsigned char)((v==255*4)?255:0);
        srcMask += 2;
      }
      // one possible unpaired column
      if (xs&1)
      {
        // image
        v = ((int)(srcDif[   0])+(int)(srcDif[xs  ])+1)/2;
        *dstPyr++ = (unsigned char)v;
        srcDif += 1;
        // mask
        v = (int)(srcMask[   0])+(int)(srcDif[xs  ]);
//        *dstMask++ = (unsigned char)((v)?255:0);
*dstMask++ = (unsigned char)((v==255*2)?255:0);
        srcMask += 1;
      }
      srcDif += xs;
      srcMask += xs;
    }
    // one possible unpaired string
    if (ys&1)
    {
      // pairs of columns
      for (j=0;j<xs/2;j++)
      {
        // image
        v = ((int)(srcDif[   0])+(int)(srcDif[   1])+1)/2;
        *dstPyr++ = (unsigned char)v;
        srcDif += 2;
        //mask
        v = (int)(srcMask[   0])+(int)(srcMask[   1]);
//        *dstMask++ = (unsigned char)((v)?255:0);
*dstMask++ = (unsigned char)((v==255*2)?255:0);
        srcMask += 2;
      }
      // one possible unpaired column
      if (xs&1)
      {
        // image
        *dstPyr++ = srcDif[   0];
        srcDif += 1;
        // mask
        *dstMask++ = srcMask[   0];
        srcMask += 1;
      }
    }
    // second pass - calculating difference histogram
    dstPyr = dstPyr_stor;
    srcDif = srcDif_stor;
    dstMask = dstMask_stor;
    srcMask = srcMask_stor;
    xsp = (xs+1)/2;
    ysp = (ys+1)/2;
    // pairs of strings
    for (i=0;i<ys/2;i++)
    {
      // pairs of columns
      for (j=0;j<xs/2;j++)
      {
        if (srcMask[0])
        {
          dx = (j       )?(   -1):0;
          dy = (i       )?(-xsp):0;
          v = (5*dstPyr[0]+dstPyr[dx+dy]+dstPyr[dx]+dstPyr[dy]+4)/8;
          buf_of_512ints[256+(v-srcDif[   0])]++;
        }
        if (srcMask[1])
        {
          dx = (j<xsp-1)?(    1):0;
          dy = (i       )?(-xsp):0;
          v = (5*dstPyr[0]+dstPyr[dx+dy]+dstPyr[dx]+dstPyr[dy]+4)/8;
          buf_of_512ints[256+(v-srcDif[   1])]++;
        }
        if (srcMask[xs])
        {
          dx = (j       )?(   -1):0;
          dy = (i<ysp-1)?( xsp):0;
          v = (5*dstPyr[0]+dstPyr[dx+dy]+dstPyr[dx]+dstPyr[dy]+4)/8;
          buf_of_512ints[256+(v-srcDif[xs  ])]++;
        }
        if (srcMask[xs+1])
        {
          dx = (j<xsp-1)?(    1):0;
          dy = (i<ysp-1)?( xsp):0;
          v = (5*dstPyr[0]+dstPyr[dx+dy]+dstPyr[dx]+dstPyr[dy]+4)/8;
          buf_of_512ints[256+(v-srcDif[xs+1])]++;
        }
        dstPyr++;
        srcDif += 2;
        dstMask++;
        srcMask += 2;
      }
      // one possible unpaired column
      if (xs&1)
      {
        if (srcMask[0])
        {
          dy = (i       )?(-xsp):0;
          v = (3*dstPyr[0]+dstPyr[dy]+2)/4;
          buf_of_512ints[256+(v-srcDif[   0])]++;
        }
        if (srcMask[xs])
        {
          dy = (i<ysp-1)?( xsp):0;
          v = (3*dstPyr[0]+dstPyr[dy]+2)/4;
          buf_of_512ints[256+(v-srcDif[xs  ])]++;
        }
        dstPyr++;
        srcDif += 1;
        dstMask++;
        srcMask += 1;
      }
      srcDif += xs;
      srcMask += xs;
    }
    // one possible unpaired string
    if (ys&1)
    {
      // pairs of columns
      for (j=0;j<xs/2;j++)
      {
        if (srcMask[0])
        {
          dx = (j       )?(   -1):0;
          v = (3*dstPyr[0]+dstPyr[dx]+2)/4;
          buf_of_512ints[256+(v-srcDif[   0])]++;
        }
        if (srcMask[1])
        {
          dx = (j<xsp-1)?(    1):0;
          v = (3*dstPyr[0]+dstPyr[dx]+2)/4;
          buf_of_512ints[256+(v-srcDif[   1])]++;
        }
        dstPyr++;
        srcDif += 2;
        dstMask++;
        srcMask += 2;
      }
    }
    // calculate difference energy
    for (i=v=0;i<512;i++)
    {
      v += buf_of_512ints[i]*(i-256)*(i-256);
      Mass += buf_of_512ints[i];
    }
    if (!Mass)
      return 0;
    return (int)(sqrt(((double)v)/Mass)*256);
}

// a = [6/((n-1)*n*(n+1))] * [2*Sxy-(n+1)*Sy]
staticFunc int MIA_MSD_CalculateAngle(
  const int* vals,      // IN:  values of function from integer argument
  int n)                // IN:  number of values
{                       // RET: <angle value>*256
  int Sy,Sxy,i;

    if (n<2)
      return 0;
    for (Sy=Sxy=i=0;i<n;i++)
    {
      Sy  += vals[i];
      Sxy += vals[i]*i;
    }
    return 100*6*(2*Sxy-(n-1)*Sy)/((n-1)*n*(n+1));
}

// estimate image sharpness by pyramidal fractal dimension
MIA_RESULT_CODE IPL_PYR_EstimateSharpness(
  int* pnSharpness,         // OUT: motion estimation
  const unsigned char* im,  // IN:  source image
  int xs,                   // IN:  source image size
  int ys,                   //
  void* pvTmpBuf,           // IN:  external buffer
  int* pnBufSiz)            // IN/OUT: bytes alloced/used in buffer
{
  int j,sz,*buf_of_512ints,xs_cur,ys_cur,*energies;
  unsigned char *procdst;
  const unsigned char *procsrc;

    if (pnBufSiz==NULL)
      return ERR_GEN_NULLPOINTER;
    if ((xs<32)||(ys<32))
      return ERR_GEN_INVALID_PARAMETER;
    // calculate size
    sz = xs*ys/2+sizeof(int)*(512+32);
    if (pvTmpBuf==NULL)
    {
      *pnBufSiz = sz;
      return ERR_OK;
    }
    if (*pnBufSiz<sz)
    {
      *pnBufSiz = sz;
      return ERR_GEN_INSUFFICIENT_BUFFER;
    }
    // assign pointers
    buf_of_512ints = (int*)pvTmpBuf;
    energies = buf_of_512ints+512;
    procdst = (unsigned char*)(energies+32);
    procsrc = im;
    // iterate
    xs_cur = xs;
    ys_cur = ys;
    for (j=0;(xs_cur>1)&&(ys_cur>1);j++)
    {
      // calculate energy and obtain next pyr level
      energies[j] = IPL_PYR_GetEnergyBLI2(
        procdst,    // OUT: downscale2-pyramid
        procsrc,    // IN/OUT: source image/difference
        xs_cur,                   // IN: source image size
        ys_cur,                   //
        buf_of_512ints);
      // normalize energy to log scale
      energies[j] = (int)(energies[j]/pow(2,((double)j)/2));
      // scale image sizes
      xs_cur = (xs_cur+1)/2;
      ys_cur = (ys_cur+1)/2;
      // take next pyr level as a source
      procsrc = procdst;
      // point destination to free space in external buffer
      procdst += xs_cur*ys_cur;
    }
    *pnSharpness = -MIA_MSD_CalculateAngle(energies,j-1);
    return ERR_OK;
}

#include <stdio.h>

// estimate image sharpness by pyramidal fractal dimension
MIA_RESULT_CODE IPL_PYR_EstimateSharpness_withMask(
  int* pnSharpness,           // OUT: motion estimation
  const unsigned char* im,    // IN:  source image
  const unsigned char* mask,  // IN:  mask
  int xs,                     // IN:  source image size
  int ys,                     //
  void* pvTmpBuf,             // IN:  external buffer
  int* pnBufSiz)              // IN/OUT: bytes alloced/used in buffer
{
  int j,sz,*buf_of_512ints,xs_cur,ys_cur,*energies;
  unsigned char *procdst,*procdstmask;
  const unsigned char *procsrc,*procsrcmask;
  int sx,n;

    if (pnBufSiz==NULL)
      return ERR_GEN_NULLPOINTER;
    if ((xs<128)||(ys<128))
      return ERR_GEN_INVALID_PARAMETER;
    // calculate size
    sz = xs*ys+sizeof(int)*(512+32);
    if (pvTmpBuf==NULL)
    {
      *pnBufSiz = sz;
      return ERR_OK;
    }
    if (*pnBufSiz<sz)
    {
      *pnBufSiz = sz;
      return ERR_GEN_INSUFFICIENT_BUFFER;
    }
    if ((im==0)||(mask==NULL)||(pnSharpness==NULL))
      return ERR_GEN_NULLPOINTER;
    // assign pointers
    buf_of_512ints = (int*)pvTmpBuf;
    energies = buf_of_512ints+512;
    procdst = (unsigned char*)(energies+32);
    procdstmask = procdst+xs*ys/2;
    procsrc = im;
    procsrcmask = mask;
    // iterate
    xs_cur = xs;
    ys_cur = ys;
    for (j=0;(xs_cur>16)&&(ys_cur>16);)
    {
/*{
  char nama[FILENAME_MAX];
  sprintf(nama,"c:\\im%02d.bmp",j);
  DBGL_FILE_SaveUint8Image(procsrc,xs_cur,xs_cur,ys_cur,nama);
  sprintf(nama,"c:\\mask%02d.bmp",j);
  DBGL_FILE_SaveUint8Image(procsrcmask,xs_cur,xs_cur,ys_cur,nama);
}*/
      // calculate energy and obtain next pyr level
      energies[j] = IPL_PYR_GetEnergyBLI2_withMask(
        procdst,    // OUT: downscale2-pyramid
        procsrc,    // IN/OUT: source image/difference
        procdstmask,
        procsrcmask,
        xs_cur,                   // IN: source image size
        ys_cur,                   //
        buf_of_512ints);
      // normalize energy to log scale
      energies[j] = (int)(energies[j]/pow(2,((double)j)/2));
      j++;
      // scale image sizes
      xs_cur = (xs_cur+1)/2;
      ys_cur = (ys_cur+1)/2;
      // take next pyr level as a source
      procsrc = procdst;
      procsrcmask = procdstmask;
      // point destination to free space in external buffer
      procdst += xs_cur*ys_cur;
      procdstmask += xs_cur*ys_cur;
    }
    {
      for (sx=n=0;n<j;n++)
        sx += energies[n];
      if (sx<=0)
        return ERR_GEN_NO_DATA;
    }
//    *pnSharpness = -MIA_MSD_CalculateAngle(energies,j-1);
    *pnSharpness = (energies[1]+energies[2])-(energies[j-2]+energies[j-1]);
//    *pnSharpness = 100*((energies[1]+energies[2])-(energies[j-2]+energies[j-1]))/sx;
    return ERR_OK;
}

staticFunc void ownPyr_SmallSobel(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys)
{
  int i,j,vx,vy;

    src += xs+1;
    for (i=1;i<ys-1;i++)
    {
      for (j=1;j<xs-1;j++)
      {
        vx = src[ xs+1]+2*src[  1]+src[-xs+1]
            -src[ xs-1]-2*src[ -1]-src[-xs-1];
        vy = src[ xs+1]+2*src[ xs]+src[ xs-1]
            -src[-xs+1]-2*src[-xs]-src[-xs-1];
        if (vx<0)
          vx = -vx;
        if (vy<0)
          vy = -vy;
        vx = (vx+vy)/8;
        if (vx>255)
          vx = 255;
        *dst++ = (unsigned char)vx;
        src++;
      }
      src += 2;
      dst += 2;
    }
}

staticFunc int ownPyr_CalcDiff(
  const unsigned char* im1,
  const unsigned char* im2,
  int xs,
  int ys)
{
  int i,j,v,s;

    s = 0;
    for (i=1;i<ys-1;i++)
    {
      for (j=1;j<xs-1;j++)
      {
        v = im1[0]-im2[0];
        if (v<0)
          v = -v;
        s += v;
        im1++;
        im2++;
      }
      im1 += 2;
      im2 += 2;
    }
    return (64*s)/((xs-2)*(ys-2));
}

// estimate image motion as difference of half-frames
MIA_RESULT_CODE IPL_PYR_EstimateMotion(
  int* pnMotion,            // OUT: motion estimation
  const unsigned char* im,  // IN:  source image
  int xs,                   // IN:  source image size
  int ys,                   //
  void* pvTmpBuf,           // IN:  external buffer
  int* pnBufSiz)            // IN/OUT: bytes alloced/used in buffer
{
  int sz;
  unsigned char *split_im,*sobel_im;

    if (pnBufSiz==NULL)
      return ERR_GEN_NULLPOINTER;
    if ((xs<32)||(ys<32))
      return ERR_GEN_INVALID_PARAMETER;
    // calculate size
    sz = xs*ys*2;
    if (pvTmpBuf==NULL)
    {
      *pnBufSiz = sz;
      return ERR_OK;
    }
    if (*pnBufSiz<sz)
    {
      *pnBufSiz = sz;
      return ERR_GEN_INSUFFICIENT_BUFFER;
    }
    // assign pointers
    split_im = (unsigned char*)pvTmpBuf;
    sobel_im = split_im+xs*ys;
    // split image
    IPL_PMP_SplitHalfFrames(split_im,im,xs,ys,0);
    // calculate Sobel gradient
    ownPyr_SmallSobel(sobel_im,split_im,xs,ys);
    // calculate difference
    *pnMotion = ownPyr_CalcDiff(sobel_im,sobel_im+(xs*ys)/2,xs,ys/2);
    return ERR_OK;
}

double ownIPL_CalcHistEpsDiscr(
  double lambda,
  const void* context)  // [0]:pointer to hist, [1]:size
{
  int i;
  double C,S,V;
  const ptr_t* __phist__ = ((ptr_t**)context)[0];
  int __sz__ = (int)(((ptr_t*)context)[1]);

    if (__phist__==NULL)
      return 1e+100;
    if (lambda<0)
      return 1e+100;
    if (lambda==0)
      C = 1./256;
    else
      C = (1.-exp(-lambda))/(1-exp(-lambda*256));
    C *= __sz__;
    S = 0;
    for (i=0;i<256;i++)
    {
      V = C*exp(-lambda*i)-__phist__[i];
      S += V*V;
    }
    return S;
}

// estimate image noise as a difference between original and its median
MIA_RESULT_CODE IPL_PYR_EstimateNoise(
  int* pnNoise,             // OUT: noise estimation
  const unsigned char* im,  // IN:  source image
  unsigned int xs,                   // IN:  source image size
  unsigned int ys,                   //
  void* pvTmpBuf,           // IN:  external buffer
  int* pnBufSiz)            // IN/OUT: bytes alloced/used in buffer
{
  ptr_t *hist;
  int v,sz;
  unsigned char *median_im,*deinterlace_im;
  //int64 Sx,Sxx;
  //double sigma;
  ptr_t anContextPlace[2];

    if (pnBufSiz==NULL)
      return ERR_GEN_NULLPOINTER;
    if ((xs<32)||(ys<32))
      return ERR_GEN_INVALID_PARAMETER;
    // calculate size
    sz = 2*xs*ys+256*sizeof(hist[0]);
    if (pvTmpBuf==NULL)
    {
      *pnBufSiz = sz;
      return ERR_OK;
    }
    if (*pnBufSiz<sz)
    {
      *pnBufSiz = sz;
      return ERR_GEN_INSUFFICIENT_BUFFER;
    }
    // assign pointers
    median_im = (unsigned char*)pvTmpBuf;
    deinterlace_im = median_im+xs*ys;
    hist = (ptr_t*)(deinterlace_im+xs*ys);
    // deinterlace
    IPL_PMP_SplitHalfFrames(deinterlace_im,im,xs,ys,0);
    // median filter
    IPL_FILT_Median3x3(median_im,deinterlace_im,xs,ys);
    // calculate difference histogram
    memset(hist,0,256*sizeof(hist[0]));
    for (sz=xs*ys-1;sz>=0;sz--)
	{
      if ((v = median_im[sz]-deinterlace_im[sz])>0)
        hist[ v]++;
      else
        hist[-v]++;
	}
//py printf("\n");
//	for(sz=0;sz<256;sz++) printf("hist[%d]=%d\n",sz,hist[sz]);
    // calculate difference
//    *pnNoise =   IPL_HIST_GetQuantile((unsigned int*)hist,xs*ys,950)-
  //             2*IPL_HIST_GetQuantile((unsigned int*)hist,xs*ys,900);

/*  Sx = Sxx = 0;
    for (v=1;v<255;v++)
    {
      Sx  += v*hist[v];
      Sxx += v*v*hist[v];
    }
    sigma = ((double)Sx)/(xs*ys);
    sigma = ((double)Sxx)/(xs*ys)-sigma*sigma;
    *pnNoise = (int)(sqrt(sigma)*1000); */
    {
      double a,b,c,eps,xmin,fmin;
      //__phist__ = &(hist[0]);
      //__sz__ = xs*ys;
      anContextPlace[0] = (ptr_t)(&(hist[0]));
      anContextPlace[1] = xs*ys;
      a = 0.;
      b = 1.;
      c = 10.;
      eps = 0.000001;
      xmin = ownIPL_CalcHistEpsDiscr(a,&(anContextPlace[0]));
      xmin = ownIPL_CalcHistEpsDiscr(b,&(anContextPlace[0]));
      xmin = ownIPL_CalcHistEpsDiscr(c,&(anContextPlace[0]));
      xmin = fmin = -1;
      sz = BPL_OPT_LocateMinimum(
        &ownIPL_CalcHistEpsDiscr,
        &(anContextPlace[0]),
        &a,
        &b,
        &c,
        &eps,
        &xmin,
        &fmin);
      //__phist__ = NULL;
      fmin = sqrt(fmin)/(xs*ys);
      *pnNoise = (int)(xmin*1000);
    }
    return ERR_OK;
}

// detect loss of sync as a spike of differences between adjacent lines
MIA_RESULT_CODE IPL_PYR_DetectSyncLoss(
  int* pnSyncLoss,          // OUT: sync loss detection flag
  const unsigned char* im,  // IN:  source image
  int xs,                   // IN:  source image size
  int ys,                   //
  void* pvTmpBuf,           // IN:  external buffer
  int* pnBufSiz)            // IN/OUT: bytes alloced/used in buffer
{
  int i,j,*hist,v,s,max0,max1,max2;

    if (pnBufSiz==NULL)
      return ERR_GEN_NULLPOINTER;
    if ((xs<32)||(ys<32))
      return ERR_GEN_INVALID_PARAMETER;
    // calculate size
    i = 256*sizeof(hist[0]);
    if (pvTmpBuf==NULL)
    {
      *pnBufSiz = i;
      return ERR_OK;
    }
    if (*pnBufSiz<i)
    {
      *pnBufSiz = i;
      return ERR_GEN_INSUFFICIENT_BUFFER;
    }
    // assign pointers
    hist = (int*)pvTmpBuf;
    // calculate difference histogram
    memset(hist,0,256*sizeof(hist[0]));
    for (i=0;i<ys-2;i+=2)
    {
      s = 0;
      for (j=0;j<xs;j++)
      {
        v = im[0]-im[2*xs];
        s += v*v;
        im++;
      }
      im += xs;
      hist[(int)sqrt(s/xs)]++;
    }
    // find maximum difference
    max0 = max1 = max2 = 255;
    for (i=255;i>=0;i--)
      if (hist[i])
      {
        max0 = i;
        hist[i]--;
        break;
      }
    // find first but maximum difference
    for (;i>=0;i--)
      if (hist[i])
      {
        max1 = i;
        hist[i]--;
        break;
      }
    // find second but maximum difference
    for (;i>=0;i--)
      if (hist[i])
      {
        max2 = i;
        break;
      }
    *pnSyncLoss = (((max0-max1)>(max1-max2))&&(max0-max1>1));
    return ERR_OK;
}
