/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  30 March 2001                                          */
/*      Revision: 1.0.00                                                 */
/*      Purpose:  operator pumps                                         */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/
#define __FILENUM__ 112 // __FILENUM__TAG112

#include <memory.h>
#include <math.h>
#include <malloc.h>
#include "ownipl.h"

// calculate X mirroring coordinate LUT
staticFunc void ownIPL_CalculateXLUT(
  int* xlut,   // OUT: X LUT
  int xs,      // IN:  width of image
  int x,       // IN:  start position of ROI
  int w,       // IN:  width of ROI
  int kew)     // IN:  kernel half-width
{
  int i;

    for (i=x-kew;i<x+w+kew;i++)
      if (i>=0)
        if (i<xs)
          xlut[i] = i;
        else
          xlut[i] = 2*xs-i-1;
      else
        xlut[i] = -(i+1);
}

// calculate Y mirroring coordinate LUT
staticFunc void ownIPL_CalculateYLUT(
  const unsigned char** ylut,  // OUT: Y LUT
  const unsigned char* src,    // IN:  source pointer
  int srcstr,                  // IN:  source stride
  int ys,                      // IN:  image height
  int y,                       // IN:  start of ROI
  int h,                       // IN: height of ROI
  int keh)                     // kernel half-height
{
  int i;

    for (i=y-keh;i<y+h+keh;i++)
      if (i>=0)
        if (i<ys)
          ylut[i] = src+srcstr*i;
        else
          ylut[i] = src+srcstr*(2*ys-i-1);
      else
        ylut[i] = src+srcstr*(-(i+1));
}

// calculate Y mirroring coordinate LUT
staticFunc void ownIPL_CalculateYLUT_int(
  const int** ylut,  // OUT: Y LUT
  const int* src,    // IN:  source pointer
  int srcstr,                  // IN:  source stride
  int ys,                      // IN:  image height
  int y,                       // IN:  start of ROI
  int h,                       // IN: height of ROI
  int keh)                     // kernel half-height
{
  int i;

    for (i=y-keh;i<y+h+keh;i++)
      if (i>=0)
        if (i<ys)
          ylut[i] = src+srcstr*i;
        else
          ylut[i] = src+srcstr*(2*ys-i-1);
      else
        ylut[i] = src+srcstr*(-(i+1));
}

#define mymax(a,b) (((a) > (b)) ? (a) : (b))
#define mymin(a,b) (((a) < (b)) ? (a) : (b))

// local normalization
MIA_RESULT_CODE IPL_FILT_LocalNormalise(
  unsigned char* dst,     // dst
  int dststr,   // dest stride
  const unsigned char* src,     // src
  int srcstr,   // source stride
  int xs,       // xs
  int ys,       // ys
  int x,        // ROI
  int y,        //
  int w,        //
  int h,        //
  int kew,      // kernel half size
  int keh)      //
{
  int xlut[PROCESSED_IMAGE_MAX_WIDTH+2*GRL_FILTER_MAX_HALFWIDTH];
  const unsigned char* ylut[PROCESSED_IMAGE_MAX_HEIGHT+2*GRL_FILTER_MAX_HALFHEIGHT];
  int sums[PROCESSED_IMAGE_MAX_WIDTH];
  int i,j,xstart,xfinish;
  int *Xlut,Sum,norm;
  const unsigned char** Ylut;

    // check arguments
    if ((xs>PROCESSED_IMAGE_MAX_WIDTH)||(ys>PROCESSED_IMAGE_MAX_HEIGHT))
      return ERR_GEN_INSUFFICIENT_BUFFER;
    // 
    norm = (2*keh+1)*(2*kew+1)*127;
    // zero-based LUT pointers
    Xlut = xlut+GRL_FILTER_MAX_HALFWIDTH;
    Ylut = ylut+GRL_FILTER_MAX_HALFHEIGHT;
    // calculate luts with mirroring
    ownIPL_CalculateXLUT(Xlut,xs,x,w,kew);
    ownIPL_CalculateYLUT(Ylut,src,srcstr,ys,y,h,keh);
    // make sums be 1 to avoid exeptions
    xstart = mymax(x-kew,0);
    xfinish = mymin(x+w+kew,xs);
    for (j=xstart;j<xfinish;j++)
      sums[j] = 1;
    for (j=xstart;j<xfinish;j++)
      for (i=y-keh;i<y+keh;i++)
        sums[j] += Ylut[i][j];
    //============ scroll rows
    for (i=y;i<y+h;i++)
    {
      // adjust sums
      xfinish = mymin(x+w+kew,xs);
      for (j=xstart;j<xfinish;j++)
        sums[j] += Ylut[i+keh][j];    // thru LUT
      // form the gross sum
      Sum = 0;
      for (j=x-kew;j<x+kew;j++)
        Sum += sums[Xlut[j]];
      //------------ left margin
      j = x;
      xfinish = mymin(kew,x+w);
      for (;j<xfinish;j++)
      {
        // adjust the gross sum
        Sum += sums[Xlut[j+kew]];    // thru LUT
        // calculate it!
        dst[i*dststr+j] = (unsigned char)((src[i*srcstr+j]*norm)/Sum);
        // adjust gross sum
        Sum -= sums[Xlut[j-kew]];    // thru LUT
      }
      //------------ center
      j = mymax(xfinish,x);
      xfinish = mymin(xs-kew,x+w);
      for (;j<xfinish;j++)
      {
        // adjust the gross sum
        Sum += sums[j+kew];    // no LUT
        // calculate it!
        dst[i*dststr+j] = (unsigned char)((src[i*srcstr+j]*norm)/Sum);
        // adjust gross sum
        Sum -= sums[j-kew];    // no LUT
      }
      //------------ right margin
      j = mymax(xfinish,x);
      xfinish = mymin(xs,x+w);
      for (;j<xfinish;j++)
      {
        // adjust the gross sum
        Sum += sums[Xlut[j+kew]];    // thru LUT
        // calculate it!
        dst[i*dststr+j] = (unsigned char)((src[i*srcstr+j]*norm)/Sum);
        // adjust gross sum
        Sum -= sums[Xlut[j-kew]];    // thru LUT
      }
      // adjust sums
      xfinish = mymin(x+w+kew,xs);
      for (j=xstart;j<xfinish;j++)
        sums[j] -= Ylut[i-keh][j];    // thru LUT
    }
    return ERR_OK;
}

// normalise globally (i.e. to fixed total brightness)
void IPL_FILT_GlobalNormalise(
  unsigned char* dst,         // destination
  int dststr,                 // stride
  const unsigned char* srcs,  // source (NULL for inplace)
  int srcstr,                 // stride
  int w,                      // ROI size
  int h,                      //
  unsigned char v)            // norm
{
  int i,j,val;
  const unsigned char* src;
  unsigned int Sum,mul;

    if (srcs==NULL)
    {
      src = (unsigned char*)dst;
      srcstr = dststr;
    }
    else
      src = srcs;
    // calculate total brightness
    Sum = 0;
    for (i=0;i<h;i++)
      for (j=0;j<w;j++)
        Sum += src[i*srcstr+j];
    // estimate multiplier
    if ((Sum>(unsigned int)(w*h*v-w*h/2))&&(Sum<(unsigned int)(w*h*v+w*h/2)))
    { // no need to change
      MIA_memcpy(dst,src,h*srcstr);
      return;
    }
    mul = w*h*(unsigned int)v;
    mul = ((mul<<8)+0x80)/Sum;
    // normalise
    for (i=0;i<h;i++)
      for (j=0;j<w;j++)
      {
        val = ((unsigned int)(src[i*srcstr+j])*mul)>>8;
        if (val<256)
          dst[i*dststr+j] = (unsigned char)val;
        else
          dst[i*dststr+j] = 255;
      }
}

// average
MIA_RESULT_CODE IPL_FILT_CalculateMean(
  unsigned char* dst,       // dst
  int dststr,              // dest stride
  const unsigned char* src, // src
  int srcstr,   // source stride
  int xs,       // xs
  int ys,       // ys
  int x,        // ROI
  int y,        //
  int w,        //
  int h,        //
  int kew,      // kernel half size
  int keh)      //
{
  int xlut[PROCESSED_IMAGE_MAX_WIDTH+2*GRL_FILTER_MAX_HALFWIDTH];
  const unsigned char* ylut[PROCESSED_IMAGE_MAX_HEIGHT+2*GRL_FILTER_MAX_HALFHEIGHT];
  int sums[PROCESSED_IMAGE_MAX_WIDTH];
  int i,j,xstart,xfinish;
  int *Xlut,Sum,norm;
  const unsigned char** Ylut;

    // check arguments
    if ((src==NULL)||(dst==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((kew<0)||(keh<0)||(x>=xs)||(y>=ys)||(x+w<0)||(y+h<0))
      return ERR_GEN_NO_DATA;
    if ((2*kew>=xs)||(2*keh>=ys)||
        (xs>PROCESSED_IMAGE_MAX_WIDTH)||(ys>PROCESSED_IMAGE_MAX_HEIGHT)||
        (kew>GRL_FILTER_MAX_HALFWIDTH)||(keh>GRL_FILTER_MAX_HALFHEIGHT)||
        (xs>srcstr)||(xs>dststr))
      return ERR_GEN_INVALID_PARAMETER;
    // process
    norm = (2*keh+1)*(2*kew+1);
    // zero-based LUT pointers
    Xlut = xlut+GRL_FILTER_MAX_HALFWIDTH;
    Ylut = ylut+GRL_FILTER_MAX_HALFHEIGHT;
    // calculate luts with mirroring
    ownIPL_CalculateXLUT(Xlut,xs,x,w,kew);
    ownIPL_CalculateYLUT(Ylut,src,srcstr,ys,y,h,keh);
    // make sums be 1 to avoid exeptions
    for (j=mymax(x-kew,0);j<mymin(x+w+kew,xs);j++)
      sums[j] = 1;
    xstart = mymax(x-kew,0);
    xfinish = mymin(x+w+kew,xs);
    for (j=xstart;j<xfinish;j++)
      for (i=y-keh;i<y+keh;i++)
        sums[j] += Ylut[i][j];
    //============ scroll rows
    for (i=y;i<y+h;i++)
    {
      // adjust sums
      xfinish = mymin(x+w+kew,xs);
      for (j=xstart;j<xfinish;j++)
        sums[j] += Ylut[i+keh][j];    // thru LUT
      // form the gross sum
      Sum = 0;
      for (j=x-kew;j<x+kew;j++)
        Sum += sums[Xlut[j]];
      //------------ left margin
      j = x;
      xfinish = mymin(kew,x+w);
      for (;j<xfinish;j++)
      {
        // adjust the gross sum
        Sum += sums[Xlut[j+kew]];    // thru LUT
        // calculate it!
        dst[i*dststr+j] = (unsigned char)(Sum/norm);
        // adjust gross sum
        Sum -= sums[Xlut[j-kew]];    // thru LUT
      }
      //------------ center
      j = mymax(xfinish,x);
      xfinish = mymin(xs-kew,x+w);
      for (;j<xfinish;j++)
      {
        // adjust the gross sum
        Sum += sums[j+kew];    // no LUT
        // calculate it!
        dst[i*dststr+j] = (unsigned char)(Sum/norm);
        // adjust gross sum
        Sum -= sums[j-kew];    // no LUT
      }
      //------------ right margin
      j = mymax(xfinish,x);
      xfinish = mymin(xs,x+w);
      for (;j<xfinish;j++)
      {
        // adjust the gross sum
        Sum += sums[Xlut[j+kew]];    // thru LUT
        // calculate it!
        dst[i*dststr+j] = (unsigned char)(Sum/norm);
        // adjust gross sum
        Sum -= sums[Xlut[j-kew]];    // thru LUT
      }
      // adjust sums
      xfinish = mymin(x+w+kew,xs);
      for (j=xstart;j<xfinish;j++)
        sums[j] -= Ylut[i-keh][j];    // thru LUT
    }
    return ERR_OK;
}

// calculation of mean and dispersion on image
MIA_RESULT_CODE IPL_FILT_CalculateMeanDisp(
  unsigned char* mean,       // dst
  unsigned char* disp,       // dst
  int dststr,              // dest stride
  const unsigned char* src, // src
  int srcstr,   // source stride
  int xs,       // xs
  int ys,       // ys
  int x,        // ROI
  int y,        //
  int w,        //
  int h,        //
  int kew,      // kernel half size
  int keh,      //
  int dispmul)  // dispersion multyplier
{
  int xlut[PROCESSED_IMAGE_MAX_WIDTH+2*GRL_FILTER_MAX_HALFWIDTH];
  const unsigned char* ylut[PROCESSED_IMAGE_MAX_HEIGHT+2*GRL_FILTER_MAX_HALFHEIGHT];
  int sums[PROCESSED_IMAGE_MAX_WIDTH],sums2[PROCESSED_IMAGE_MAX_WIDTH];
  int i,j,xstart,xfinish,v;
  int *Xlut,Sum,norm,Sum2;
  const unsigned char** Ylut;

    // check arguments
    if ((src==NULL)||(mean==NULL)||(disp==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((kew<=0)||(keh<=0)||(x>=xs)||(y>=ys)||(x+w<0)||(y+h<0))
      return ERR_GEN_NO_DATA;
    if ((2*kew>=xs)||(2*keh>=ys)||
        (xs>PROCESSED_IMAGE_MAX_WIDTH)||(ys>PROCESSED_IMAGE_MAX_HEIGHT)||
        (kew>GRL_FILTER_MAX_HALFWIDTH)||(keh>GRL_FILTER_MAX_HALFHEIGHT)||
        (xs>srcstr)||(xs>dststr))
      return ERR_GEN_INVALID_PARAMETER;
    // process
    norm = (2*keh+1)*(2*kew+1);
    // zero-based LUT pointers
    Xlut = xlut+GRL_FILTER_MAX_HALFWIDTH;
    Ylut = ylut+GRL_FILTER_MAX_HALFHEIGHT;
    // calculate luts with mirroring
    ownIPL_CalculateXLUT(Xlut,xs,x,w,kew);
    ownIPL_CalculateYLUT(Ylut,src,srcstr,ys,y,h,keh);
    for (j=mymax(x-kew,0);j<mymin(x+w+kew,xs);j++)
      sums[j] = sums2[j] = 0;
    xstart = mymax(x-kew,0);
    xfinish = mymin(x+w+kew,xs);
    for (j=xstart;j<xfinish;j++)
      for (i=y-keh;i<y+keh;i++)
      {
        v = Ylut[i][j];
        sums[j] += v;
        sums2[j] += v*v;
      }
    //============ scroll rows
    for (i=y;i<y+h;i++)
    {
      // adjust sums
      xfinish = mymin(x+w+kew,xs);
      for (j=xstart;j<xfinish;j++)
      {
        v = Ylut[i+keh][j]; // thru LUT
        sums[j] += v;
        sums2[j] += v*v;
      }
      // form the gross sum
      Sum = Sum2 = 0;
      for (j=x-kew;j<x+kew;j++)
      {
        Sum += sums[Xlut[j]];
        Sum2 += sums2[Xlut[j]];
      }
      //------------ left margin
      j = x;
      xfinish = mymin(kew,x+w);
      for (;j<xfinish;j++)
      {
        // adjust the gross sum
        Sum += sums[Xlut[j+kew]];    // thru LUT
        Sum2 += sums2[Xlut[j+kew]];    // thru LUT
        // calculate it!
        mean[i*dststr+j] = (unsigned char)(Sum/norm);
        v = (int)(dispmul*sqrt(((double)Sum2-((double)Sum*Sum)/norm)/norm)+.5);
        if (v>255)
          v = 255;
        disp[i*dststr+j] = (unsigned char)v;
        // adjust gross sum
        Sum -= sums[Xlut[j-kew]];    // thru LUT
        Sum2 -= sums2[Xlut[j-kew]];    // thru LUT
      }
      //------------ center
      j = mymax(xfinish,x);
      xfinish = mymin(xs-kew,x+w);
      for (;j<xfinish;j++)
      {
        // adjust the gross sum
        Sum += sums[j+kew];    // no LUT
        Sum2 += sums2[j+kew];    // no LUT
        // calculate it!
        mean[i*dststr+j] = (unsigned char)(Sum/norm);
        v = (int)(dispmul*sqrt(((double)Sum2-((double)Sum*Sum)/norm)/norm)+.5);
        if (v>255)
          v = 255;
        disp[i*dststr+j] = (unsigned char)v;
        // adjust gross sum
        Sum -= sums[j-kew];    // no LUT
        Sum2 -= sums2[j-kew];    // no LUT
      }
      //------------ right margin
      j = mymax(xfinish,x);
      xfinish = mymin(xs,x+w);
      for (;j<xfinish;j++)
      {
        // adjust the gross sum
        Sum += sums[Xlut[j+kew]];    // thru LUT
        Sum2 += sums2[Xlut[j+kew]];    // thru LUT
        // calculate it!
        mean[i*dststr+j] = (unsigned char)(Sum/norm);
        v = (int)(dispmul*sqrt(((double)Sum2-((double)Sum*Sum)/norm)/norm)+.5);
        if (v>255)
          v = 255;
        disp[i*dststr+j] = (unsigned char)v;
        // adjust gross sum
        Sum -= sums[Xlut[j-kew]];    // thru LUT
        Sum2 -= sums2[Xlut[j-kew]];    // thru LUT
      }
      // adjust sums
      xfinish = mymin(x+w+kew,xs);
      for (j=xstart;j<xfinish;j++)
      {
        v = Ylut[i-keh][j];    // thru LUT
        sums[j] -= v;
        sums2[j] -= v*v;
      }
    }
    return ERR_OK;
}

float IPL_FILT_Mean_Disp(
  short *bufin,
  short *Midl,
  int lstr)
{
  int i,ts,Num=0;
  int S=0,M=0;
  float D;

    for (i=0;i<lstr;) 
      if((ts=bufin[i++])>0)
      {
        M+=ts;
        S+=ts*ts;
        Num++;
      }
    if(Num>0)
    {
      *Midl=(short)((M+Num/2)/Num);
      D=(float)(Num*S-M*M);
      D=(float)(sqrt(D)/Num);
      return D;
    }
    *Midl=0;
    return 1000;
}

// calculate local variation
MIA_RESULT_CODE IPL_FILT_LocalVariation(
  unsigned char* dst,
  const unsigned char* src,
  unsigned int xs,
  unsigned int ys,
  unsigned int wx,
  unsigned int wy)
{
  unsigned int  _hist[258];
  unsigned int * hist = _hist+1;
  unsigned int  _min_,_max_,p;
  int i,j;
  unsigned char *finish;

    memset(_hist,0,258*sizeof(int));
    // initial fill and find
    finish = dst+(ys-(wy+1))*xs+(xs-(wx+1)); // finishing pointer
    src += wy*xs+wx;  // starting at point (wx,wy)
    dst += wy*xs+wx; // finishing at (xs-(wx+1),ys-(wy+2))
    _min_ = 255;
    _max_ = 0;
    for (i=-(int)(wy);i<=(int)(wy);i++)     // init. histogram
      for (j=-(int)(wx);j<=(int)(wx);j++)
      {
        p = src[i*xs+j];
        if (p>_max_)
          _max_ = p;
        if (p<_min_)
          _min_ = p;
        hist[p]++;
      }
    j = wy*xs-wx; // shortcut
    while (finish!=dst)
    {
      *dst++ = (unsigned char)(_max_-_min_);
      for (i=-(int)(wy)*xs-(int)(wx);i<=j;i+=xs)
      {
        p = src[i+(wx<<1)+1];
        hist[p]++;
        if (p<_min_)
          _min_ = p;
        else
          if (p>_max_)
            _max_ = p;
        p = src[i];
        if (--hist[p]==0)
        {
          if (p==_min_)
            while (hist[_min_]==0)
              _min_++;
          else
            if (p==_max_)
              while (hist[_max_]==0)
                _max_--;
        }
      }
      src++;
    }
    return ERR_OK;
}

// calculate val-mymin in local nbhood
MIA_RESULT_CODE IPL_FILT_LocalMaxDrop(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  int wx,
  int wy)
{
  int  _hist[258];
  int * hist = _hist+1;
  int  _min_,p;
  int i,j;
  unsigned char *finish,*dst_sav=dst;

    memset(_hist,0,258*sizeof(int));
    // initial fill and find
    finish = dst+(ys-(wy+1))*xs+(xs-(wx+1)); // finishing pointer
    src += wy*xs+wx;  // starting at point (wx,wy)
    dst += wy*xs+wx; // finishing at (xs-(wx+1),ys-(wy+2))
    _min_ = 255;
    for (i=-wy;i<=wy;i++)     // init. histogram
      for (j=-wx;j<=wx;j++)
      {
        p = src[i*xs+j];
        if (p<_min_)
          _min_ = p;
        hist[p]++;
      }
    j = wy*xs-wx; // shortcut
    while (finish!=dst)
    {
      // BUGBUG in case src=255, mymin=0:   0 is outputted instead of 256.
      //        Still, can be neglected in real images.
      *dst++ = (unsigned char)((*src>=_min_)?(1L+*src-_min_):0);
      for (i=-wy*xs-wx;i<=j;i+=xs)
      {
        p = src[i+(wx<<1)+1];
        hist[p]++;
        if (p<_min_)
          _min_ = p;
        p = src[i];
        if (--hist[p]==0)
          if (p==_min_)
            while (hist[_min_]==0)
              _min_++;
      }
      src++;
    }
    ownIPL_FillBordersXB(dst_sav,xs,ys,wx,wy);
    return ERR_OK;
}

// finds chessboard crosses with black quadrants 1 and 3
// processes full image, outputs full image
MIA_RESULT_CODE IPL_FILT_ChessCrossesWithBlack13Quadrant(
  unsigned char* dst,         // destination image, stride=xs
  const unsigned char* er,    // source image (eroded)
  const unsigned char* di,    // source image (dilated)
  unsigned int srcstr,       // its stride
  unsigned int xs,           // image width
  unsigned int ys,           // image height
  unsigned int kew)          // kernel half-size
{
  unsigned int add2src,i,dy;
  int val_m,val_t,val_M;

    // check arguments
    if ((dst==NULL)||(er==NULL)||(di==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((xs==0)||(ys==0)||(srcstr==0))
      return ERR_GEN_NO_DATA;
    if (kew==0)
      return ERR_GEN_INVALID_PARAMETER;
    // clear borders
    memset(dst,0,kew*xs);
    memset(dst+(ys-kew)*xs,0,kew*xs);
    dst += kew*(xs-1);
    for (i=ys-2*kew+1;i;i--)
    {
      memset(dst,0,kew*2);
      dst += xs;
    }
    // put pointers to ROI start
    dy = srcstr*(kew *= 2);
    dst -= (ys -= kew)*xs+(xs-kew);
    add2src = srcstr-(xs -= kew);
    // go!
    while (ys)
    {
      for (i=xs;i;i--)
      {
// a, b, c, and d are values in nodes of kernel rectangle:
//
//    a---b
//    |   |
//    c---d
//
// D = mymin{(am-bM),(am-cM),(dm-bM),(dm-cM)}-mymax{|am-dM|,|bm-cM|,|aM-dm|,|bM-cm|}
// destination = mymax{D,0}
// m <-> eroded, M <-> dilated

        // mymin{...}
        val_m            = (int)(er[0     ])-(int)(di[kew]);    // a-b
        if (val_m>(val_t = (int)(er[0     ])-(int)(di[dy ])))   // a-c
          val_m = val_t;
        if (val_m>(val_t = (int)(er[dy+kew])-(int)(di[dy ])))   // d-c
          val_m = val_t;
        if (val_m>(val_t = (int)(er[dy+kew])-(int)(di[kew])))   // d-b
          val_m = val_t;
        // mymax{...}
        val_M            = (int)(er[dy    ])-(int)(di[kew] );   // cm-bM
        if (val_M<0)
          val_M = -val_M;
        val_t            = (int)(di[dy    ])-(int)(er[kew] );   // cM-bm
        if (val_t<0)
          val_t = -val_t;
        if (val_M<val_t)
          val_M = val_t;
        val_t            = (int)(er[dy+kew])-(int)(di[0]  );    // dm-aM
        if (val_t<0)
          val_t = -val_t;
        if (val_M<val_t)
          val_M = val_t;
        val_t            = (int)(di[dy+kew])-(int)(er[0]  );    // dM-am
        if (val_t<0)
          val_t = -val_t;
        if (val_M<val_t)
          val_M = val_t;
        // mymin-mymax
        if ((val_t = val_m-val_M)<=0)
          *dst++ = 0;
        else
          *dst++ = (unsigned char)val_t;
        er++;
        di++;
      }
      er += add2src;
      di += add2src;
      dst += kew;
      ys--;
    }
    return ERR_OK;
}

// finds chessboard crosses with black quadrants 2 and 4
// processes full image, outputs full image
MIA_RESULT_CODE IPL_FILT_ChessCrossesWithBlack24Quadrant(
  unsigned char* dst,         // destination image, stride=xs
  const unsigned char* er,    // source image (eroded)
  const unsigned char* di,    // source image (dilated)
  unsigned int srcstr,       // its stride
  unsigned int xs,           // image width
  unsigned int ys,           // image height
  unsigned int kew)          // kernel half-size
{
  unsigned int add2src,i,dy;
  int val_m,val_t,val_M;

    // check arguments
    if ((dst==NULL)||(er==NULL)||(di==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((xs==0)||(ys==0)||(srcstr==0))
      return ERR_GEN_NO_DATA;
    if (kew==0)
      return ERR_GEN_INVALID_PARAMETER;
    // clear borders
    memset(dst,0,kew*xs);
    memset(dst+(ys-kew)*xs,0,kew*xs);
    dst += kew*(xs-1);
    for (i=ys-2*kew+1;i;i--)
    {
      memset(dst,0,kew*2);
      dst += xs;
    }
    // put pointers to ROI start
    dy = srcstr*(kew *= 2);
    dst -= (ys -= kew)*xs+(xs-kew);
    add2src = srcstr-(xs -= kew);
    // go!
    while (ys)
    {
      for (i=xs;i;i--)
      {
        // mymin{...}
        val_m            = (int)(er[kew])-(int)(di[0     ]);    // a-b
        if (val_m>(val_t = (int)(er[dy ])-(int)(di[0     ])))   // a-c
          val_m = val_t;
        if (val_m>(val_t = (int)(er[dy ])-(int)(di[dy+kew])))   // d-c
          val_m = val_t;
        if (val_m>(val_t = (int)(er[kew])-(int)(di[dy+kew])))   // d-b
          val_m = val_t;
        // mymax{...}
        val_M            = (int)(er[kew])-(int)(di[dy    ] );   // cm-bM
        if (val_M<0)
          val_M = -val_M;
        val_t            = (int)(di[kew])-(int)(er[dy    ] );   // cM-bm
        if (val_t<0)
          val_t = -val_t;
        if (val_M<val_t)
          val_M = val_t;
        val_t            = (int)(er[0])-(int)(di[dy+kew]  );    // dm-aM
        if (val_t<0)
          val_t = -val_t;
        if (val_M<val_t)
          val_M = val_t;
        val_t            = (int)(di[0])-(int)(er[dy+kew]  );    // dM-am
        if (val_t<0)
          val_t = -val_t;
        if (val_M<val_t)
          val_M = val_t;
        // mymin-mymax
        if ((val_t = val_m-val_M)<=0)
          *dst++ = 0;
        else
          *dst++ = (unsigned char)val_t;
        er++;
        di++;
      }
      er += add2src;
      di += add2src;
      dst += kew;
      ys--;
    }
    return ERR_OK;
}

void IPL_FILT_HaussBlur3x3(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys)
{
  int i,j;

    for (i=1;i<ys-1;i++)
    {
      for (j=1;j<xs-1;j++)
        dst[i*xs+j] = (unsigned char)(
                      ( 1*src[(i-1)*xs+(j-1)]+2*src[(i-1)*xs+(j  )]+1*src[(i-1)*xs+(j+1)]+
                        2*src[(i  )*xs+(j-1)]+4*src[(i  )*xs+(j  )]+2*src[(i  )*xs+(j+1)]+
                        1*src[(i+1)*xs+(j-1)]+2*src[(i+1)*xs+(j  )]+1*src[(i+1)*xs+(j+1)] )/16);
    }
    ownIPL_FillBorders1B(dst,xs,ys);
}

void IPL_FILT_HaussBlur3x3_Int(
  int* dst,
  const int* src,
  int xs,
  int ys,
  int mode) // 0 - normalized, 1- no normalization
{
  int i;

    if (mode)
    { // no normalization
      // body
      dst += xs+1;
      src += xs+1;
      for (i=(ys-2)*xs-2;i>0;i--)
      {
        *dst = 1*src[-xs-1]+2*src[-xs]+1*src[-xs+1]+
               2*src[   -1]+4*src[  0]+2*src[    1]+
               1*src[ xs-1]+2*src[ xs]+1*src[ xs+1];
        src++;
        dst++;
      }
      dst -= xs*ys-xs-1;
    }
    else
    { // normalized
      // body
      dst += xs+1;
      src += xs+1;
      for (i=(ys-2)*xs-2;i>0;i--)
      {
        *dst = ( 1*src[-xs-1]+2*src[-xs]+1*src[-xs+1]+
                 2*src[   -1]+4*src[  0]+2*src[    1]+
                 1*src[ xs-1]+2*src[ xs]+1*src[ xs+1] )/16;
        src++;
        dst++;
      }
      dst -= xs*ys-xs-1;
    }
    // borders
    ownIPL_FillBorders1B_int(dst,xs,ys);
}

void IPL_FILT_UniformBlur3x3_Int(
  int* dst,
  const int* src,
  int xs,
  int ys,
  int mode) // 0 - normalized, 1- no normalization
{
  int i;

    if (mode)
    { // no normalization
      dst += xs+1;
      src += xs+1;
      for (i=(ys-2)*xs-2;i>0;i--)
      {
        *dst = 1*src[-xs-1]+1*src[-xs]+1*src[-xs+1]+
               1*src[   -1]+1*src[  0]+1*src[    1]+
               1*src[ xs-1]+1*src[ xs]+1*src[ xs+1];
        src++;
        dst++;
      }
      dst -= xs*ys-xs-1;
    }
    else
    { // normalized
      dst += xs+1;
      src += xs+1;
      for (i=(ys-2)*xs-2;i>0;i--)
      {
        *dst = ( 1*src[-xs-1]+1*src[-xs]+1*src[-xs+1]+
                 1*src[   -1]+1*src[  0]+1*src[    1]+
                 1*src[ xs-1]+1*src[ xs]+1*src[ xs+1] )/9;
        src++;
        dst++;
      }
      dst -= xs*ys-xs-1;
    }
    // borders
    ownIPL_FillBorders1B_int(dst,xs,ys);
}

void IPL_FILT_HaussBlur5x5(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys)
{
  int i,j;

    for (i=2;i<ys-2;i++)
    {
      for (j=2;j<xs-2;j++)
        dst[i*xs+j] = (unsigned char)(
                      (   2*src[(i-2)*xs+(j-2)]+  7*src[(i-2)*xs+(j-1)]+ 12*src[(i-2)*xs+(j  )]+  7*src[(i-2)*xs+(j+1)]+  2*src[(i-2)*xs+(j+2)]+
                          7*src[(i-1)*xs+(j-2)]+ 31*src[(i-1)*xs+(j-1)]+ 52*src[(i-1)*xs+(j  )]+ 31*src[(i-1)*xs+(j+1)]+  7*src[(i-1)*xs+(j+2)]+
                         12*src[(i  )*xs+(j-2)]+ 52*src[(i  )*xs+(j-1)]+127*src[(i  )*xs+(j  )]+ 52*src[(i  )*xs+(j+1)]+ 12*src[(i  )*xs+(j+2)]+
                          7*src[(i+1)*xs+(j-2)]+ 31*src[(i+1)*xs+(j-1)]+ 52*src[(i+1)*xs+(j  )]+ 31*src[(i+1)*xs+(j+1)]+  7*src[(i+1)*xs+(j+2)]+
                          2*src[(i+2)*xs+(j-2)]+  7*src[(i+2)*xs+(j-1)]+ 12*src[(i+2)*xs+(j  )]+  7*src[(i+2)*xs+(j+1)]+  2*src[(i+2)*xs+(j+2)] )/571);
    }
    ownIPL_FillBordersXB(dst,xs,ys,2,2);
}

void IPL_FILT_HaussBlur1x5(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys)
{
  int i,j;

    for (i=2;i<ys-2;i++)
    {
      for (j=0;j<xs-0;j++)
        dst[i*xs+j] = (unsigned char)( (  12*src[(i-2)*xs+(j  )]+52*src[(i-1)*xs+(j  )]+127*src[(i  )*xs+(j  )]+52*src[(i+1)*xs+(j  )]+12*src[(i+2)*xs+(j  )] +125)/255);
        //dst[i*xs+j] = (int)( src[(i  )*xs+(j  )]);
    }
    ownIPL_FillBordersXB(dst,xs,ys,0,2);
}

void IPL_FILT_HaussBlur5x1(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys)
{
  int i,j;

    for (i=0;i<ys-0;i++)
    {
      for (j=2;j<xs-2;j++)
        dst[i*xs+j] = (unsigned char)( (  12*src[(i  )*xs+(j-2)]+52*src[(i  )*xs+(j-1)]+127*src[(i  )*xs+(j  )]+52*src[(i  )*xs+(j+1)]+12*src[(i  )*xs+(j+2)] )/255);
        //dst[i*xs+j] = (int)( src[(i  )*xs+(j  )]);
    }
    ownIPL_FillBordersXB(dst,xs,ys,2,0);
}

void IPL_FILT_HaussBlur1x3(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys)
{
  int i,j;

    for (i=1;i<ys-1;i++)
    {
      for (j=0;j<xs-0;j++)
        dst[i*xs+j] = (unsigned char)( (  52*src[(i-1)*xs+(j)]+127*src[(i  )*xs+(j  )]+52*src[(i+1)*xs+(j  )] )/231);
        //dst[i*xs+j] = (int)( src[(i  )*xs+(j  )]);
    }
    ownIPL_FillBordersXB(dst,xs,ys,0,1);
}

void IPL_FILT_HaussBlur3x1(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys)
{
  int i,j;

    for (i=0;i<ys-0;i++)
    {
      for (j=1;j<xs-1;j++)
        dst[i*xs+j] = (unsigned char)( (  52*src[(i)*xs+(j-1  )]+127*src[(i  )*xs+(j  )]+52*src[(i)*xs+(j+1  )] )/231);
        //dst[i*xs+j] = (int)( src[(i  )*xs+(j  )]);
    }
    ownIPL_FillBordersXB(dst,xs,ys,1,0);
}

void IPL_FILT_ShiftedDiff1x5(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  int shift)                  // shift for difference: simulate gabor filter
{
  int i,j;
  int shft_up = (shift-1)/2;
  int shft_dn = shift-shft_up;
  int tmprez = 0;

  for (i=shft_up;i<ys-shft_dn;i++)
  {
    for (j=0;j<xs-0;j++)
    {
      tmprez = 127+4*src[(i-shft_up)*xs+(j  )]-4*src[(i+shft_dn)*xs+(j  )];
      tmprez = mymin(250,tmprez);
      tmprez = mymax(5,tmprez);
      dst[i*xs+j] = (unsigned char)tmprez;
    }
  }
  ownIPL_FillBordersXB(dst,xs,ys,0,shft_dn);
}

void IPL_FILT_ShiftedDiff5x1(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  int shift)                  // shift for difference: simulate gabor filter
{
  int i,j;
  int shft_up = (shift-1)/2;
  int shft_dn = shift-shft_up;
  int tmprez = 0;

  for (i=0;i<ys;i++)
  {
    for (j=shft_up;j<xs-shft_dn;j++)
    {
      tmprez = 127+4*src[(i  )*xs+(j-shft_up)]-4*src[(i  )*xs+(j+shft_dn)];
      tmprez = mymin(250,tmprez);
      tmprez = mymax(5,tmprez);
      dst[i*xs+j] = (unsigned char)tmprez;
    }
  }
  ownIPL_FillBordersXB(dst,xs,ys,shft_dn,0);
}

void IPL_FILT_HaussBlur5x5_int32(
  int* dst,
  const int* src,
  int xs,
  int ys,
  int bNormalise)
{
  int i,j;

    // body
    if (bNormalise)
      for (i=2;i<ys-2;i++)
        for (j=2;j<xs-2;j++)
          dst[i*xs+j] = (   2*src[(i-2)*xs+(j-2)]+  7*src[(i-2)*xs+(j-1)]+ 12*src[(i-2)*xs+(j  )]+  7*src[(i-2)*xs+(j+1)]+  2*src[(i-2)*xs+(j+2)]+
                            7*src[(i-1)*xs+(j-2)]+ 31*src[(i-1)*xs+(j-1)]+ 52*src[(i-1)*xs+(j  )]+ 31*src[(i-1)*xs+(j+1)]+  7*src[(i-1)*xs+(j+2)]+
                           12*src[(i  )*xs+(j-2)]+ 52*src[(i  )*xs+(j-1)]+127*src[(i  )*xs+(j  )]+ 52*src[(i  )*xs+(j+1)]+ 12*src[(i  )*xs+(j+2)]+
                            7*src[(i+1)*xs+(j-2)]+ 31*src[(i+1)*xs+(j-1)]+ 52*src[(i+1)*xs+(j  )]+ 31*src[(i+1)*xs+(j+1)]+  7*src[(i+1)*xs+(j+2)]+
                            2*src[(i+2)*xs+(j-2)]+  7*src[(i+2)*xs+(j-1)]+ 12*src[(i+2)*xs+(j  )]+  7*src[(i+2)*xs+(j+1)]+  2*src[(i+2)*xs+(j+2)] )/571;
    else
      for (i=2;i<ys-2;i++)
        for (j=2;j<xs-2;j++)
          dst[i*xs+j] =     2*src[(i-2)*xs+(j-2)]+  7*src[(i-2)*xs+(j-1)]+ 12*src[(i-2)*xs+(j  )]+  7*src[(i-2)*xs+(j+1)]+  2*src[(i-2)*xs+(j+2)]+
                            7*src[(i-1)*xs+(j-2)]+ 31*src[(i-1)*xs+(j-1)]+ 52*src[(i-1)*xs+(j  )]+ 31*src[(i-1)*xs+(j+1)]+  7*src[(i-1)*xs+(j+2)]+
                           12*src[(i  )*xs+(j-2)]+ 52*src[(i  )*xs+(j-1)]+127*src[(i  )*xs+(j  )]+ 52*src[(i  )*xs+(j+1)]+ 12*src[(i  )*xs+(j+2)]+
                            7*src[(i+1)*xs+(j-2)]+ 31*src[(i+1)*xs+(j-1)]+ 52*src[(i+1)*xs+(j  )]+ 31*src[(i+1)*xs+(j+1)]+  7*src[(i+1)*xs+(j+2)]+
                            2*src[(i+2)*xs+(j-2)]+  7*src[(i+2)*xs+(j-1)]+ 12*src[(i+2)*xs+(j  )]+  7*src[(i+2)*xs+(j+1)]+  2*src[(i+2)*xs+(j+2)] ;
    ownIPL_FillBordersXB_int(dst,xs,ys,2,2);
}

void IPL_FILT_ClearAlonePoints8Conn(
  unsigned char* im,
  int xs,
  int ys)
{
  int i,j,s;

    for(i=ys-2;i>=1;i--)
      for(j=xs-2;j>=1;j--)
      {
        s = im[(i-1)*xs+(j-1)]+
            im[(i-1)*xs+(j  )]+
            im[(i-1)*xs+(j+1)]+
            im[(i  )*xs+(j-1)]+
            //im[(i  )*xs+(j  )]+
            im[(i  )*xs+(j+1)]+
            im[(i+1)*xs+(j-1)]+
            im[(i+1)*xs+(j  )]+
            im[(i+1)*xs+(j+1)];
        if (!s)
          im[(i  )*xs+(j  )] = 0;
      }
}

void IPL_FILT_Median3x1_Int(
  int* dst,
  const int* src,
  int xs,
  int ys)
{
  int i;
  register int j,m,c,p;

    for (i=0;i<ys;i++)
    {
      // first byte
      dst[0] = ((m = src[0])+(c = src[1]))/2;
      for (j=1;j<xs-1;j++)
      {
        p = src[j+1];
        if (m>p)
        {
          if (c<p)
            dst[j] = p;
          else
            if (m<c)
              dst[j] = m;
            else
              dst[j] = c;
        }
        else
          if (c>p)
            dst[j] = p;
          else
            if (m>c)
              dst[j] = m;
            else
              dst[j] = c;
        m = c;
        c = p;
      }
      // last byte
      dst[j] = (m+c)/2;
      src += xs;
      dst += xs;
    }
}
