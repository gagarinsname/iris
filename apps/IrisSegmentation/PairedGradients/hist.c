/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  4 August 2010                                          */
/*      Revision: 1.0.00                                                 */
/*      Purpose:  histogram processing functions                         */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#define __FILENUM__ 99 // __FILENUM__TAG99

#include <memory.h>
#include <malloc.h>
#include "ipl.h"

int IPL_HIST_GetQuantile(     // RET: position of quantile
  const unsigned int* pnHist, // IN:  histogram
  unsigned int nMass,         // IN:  histogam mass
  unsigned int nQuantile)     // IN:  quantile in 1/1000
{
  int nSum,i;

    if ((pnHist==NULL)||(nQuantile>1000))
      return 0;
    if (nMass<=0)
    {
      nMass = 0;
      for (i=0;i<256;i++)
        nMass += pnHist[i];
    }
    nSum = nMass*nQuantile/1000;
    for (i=0;(i<256)&&(nSum>0);)
      nSum -= pnHist[i++];
    if (nSum)
      return i-1;
    // nSum=0 - should return middle of (i-1) and next non-zero
    nSum = i;
    while ((i<256)&&(!pnHist[i++]));
    return (nSum+i-2)/2;
}

int IPL_HIST_GetQuantileExt(  // RET: position of quantile
  const unsigned int* pnHist, // IN:  histogram
  int nLen,           // IN:  histogram length
  int nMass,          // IN:  histogam mass
  int nQuantile)      // IN:  quantile in 1/1000
{
  int nSum,i;

    if ((pnHist==NULL)||(nQuantile>1000)||(nLen<=0))
      return 0;
    if (nMass<=0)
      for (i=0;i<nLen;i++)
        nMass += pnHist[i];
    nSum = nMass*nQuantile/1000;
    for (i=0;(i<nLen)&&(nSum>0);)
      nSum -= pnHist[i++];
    if (nSum)
      return i-1;
    // nSum=0 - should return middle of (i-1) and next non-zero
    nSum = i;
    while ((i<nLen)&&(!pnHist[i++]));
    return (nSum+i-2)/2;
}

void IPL_HIST_Blur(
  int* dst,         // OUT: blurred hist
  const int* val,   // IN:  source hist
  int len,          // IN:  length of the hist
  int hw)           // IN:  half-window
{
  int i,V;

    // beginning
    V = val[0]*(hw+1);
    for (i=1;i<hw;i++)
      V += val[i];
    for (i=0;i<hw;i++)
    {
      V += val[i+hw];
      dst[i] = V;
      V -= val[0];
    }
    // body
    len -= hw;
    for (;i<len;i++)
    {
      V += val[i+hw];
      dst[i] = V;
      V -= val[i-hw];
    }
    // ending
    len += hw;
    for (;i<len;i++)
    {
      V += val[len-1];
      dst[i] = V;
      V -= val[i-hw];
    }
}

void IPL_HIST_Blur_double(
  double* dst,         // OUT: blurred hist
  const double* val,   // IN:  source hist
  int len,          // IN:  length of the hist
  int hw)           // IN:  half-window
{
  int i;
  double V;

    // beginning
    V = val[0]*(hw+1);
    for (i=1;i<hw;i++)
      V += val[i];
    for (i=0;i<hw;i++)
    {
      V += val[i+hw];
      dst[i] = V;
      V -= val[0];
    }
    // body
    len -= hw;
    for (;i<len;i++)
    {
      V += val[i+hw];
      dst[i] = V;
      V -= val[i-hw];
    }
    // ending
    len += hw;
    for (;i<len;i++)
    {
      V += val[len-1];
      dst[i] = V;
      V -= val[i-hw];
    }
}

// needed for filling brightness histogram gaps
// occuring after color adjustment procedures
// uses linear interpolation by non-zero neighbors
void IPL_HIST_FillGaps(
  int* dst,     // IN/OUT: histogram (transformed in-place)
  int len,      // IN:  histogram data length
  int wid)      // IN:  maximum width of gap
{
  int i,id_l=0,j;

    // wind up to first non-zero
    for (i=0;(i<len)&&(dst[i]==0);i++);
    // till the end
    while (i+wid<len)
      if (dst[i])
      { // wind up to zero.
        for (;(i+wid<len)&&(dst[i]!=0);i++);
        // it'll be one of interpolative values for forthcoming gap
        id_l = i-1;
      }
      else
      {
        // wind up to non-zero.
        for (;(i<len)&&(dst[i]==0);i++);
        // is the gap smaller than maximum width?
        if (i-id_l<=wid+1)
          // yes. do interpolation
          for (j=id_l+1;j<i;j++)
            dst[j] = ((j-id_l)*dst[i]+(i-j)*dst[id_l])/(i-id_l);
      }
}

// make median of all array for byte
int IPL_HIST_mdn_pixE(
  const unsigned char* bufin, // input array
  int                  lstr)  // size bufin
{                             // Returns: median
  int i,half,S,anARM[256];

    memset(&anARM[0],0,sizeof(anARM));
    half = (lstr/2)|1;
    for (i=0;i<lstr;anARM[bufin[i++]]++);
    for (i=0,S=0;S<half;S+=anARM[i++]);
    return i-1;
}

// Returns: median
int IPL_HIST_grouping1(
  short *bufin,
  int med,
  int lstr)
{
  int i,i0,half,S,S0;
  short anARM[256],median;  

    MIA_memset(&anARM[0],0,sizeof(anARM));
    half = (lstr/2)|1;
    median= (short)med;
    for (i=0;i<lstr;i++)
    {
      S=bufin[i]-median;
      if(S<0)
        S=-S;
      anARM[S]+=1;
    }
    for (i0=0;i0<50;i0++)
      if(anARM[i0]!=0)
        break;
    for (i0=0,S0=0;S0<half;S0+=anARM[i0++]);
    for (i=i0,S=S0;S<lstr;S+=anARM[i++]);
    S0=half;
    S=lstr-S0;
    S=S*(i-i0);
//    S0=S0*32/i0;
//    return S*16/S0;
    return S*i0/(S0*2);
}

// Returns: median
int IPL_HIST_grouping(
  short *bufin,
  int med,
  int lstr)
{
  int i,half,S;
  short anARM[256],median;  

    MIA_memset(&anARM[0],0,sizeof(anARM));
    half = (lstr/2)|1;
    median= (short)med;
    for (i=0;i<lstr;i++)
    {
      S=bufin[i]-median;
      if(S<0)
        S=-S;
      anARM[S]+=1;
    }
    for (i=0,S=0;S<half;S+=anARM[i++]);
    return i-1;
}

MIA_RESULT_CODE IPL_HIST_CalcMedianInLine(
  unsigned char* line,
  const unsigned char* im,
  int xs,
  int ys,
  int h,
  int kew)
{
  int i,hist[256],vmed,dx,dy,dxb,nlef,nrig,v;

    // check arguments
    if (2*kew+1>ys)
      return ERR_GEN_NO_DATA;
    if (h<kew)
      h = kew;
    if (h+kew>=ys)
      h = ys-kew-1;
    // initial set
    im += h*xs;
    memset(hist,0,sizeof(hist));
    for (dx=-kew;dx<kew;dx++) // '<' rather than '<='
    {
      dxb = MIA_bound(dx,0,xs-1);
      for (dy=-kew;dy<=kew;dy++)
        hist[im[dy*xs+dxb]]++;
    }
    vmed = 0;
    nlef = 0;
    nrig = (2*kew+1)*(2*kew)-hist[0];
    // main cycle
    for (i=0;i<xs;i++)
    {
      // add points from right
      dxb = MIA_bound(i+kew,0,xs-1);
      for (dy=-kew;dy<=kew;dy++)
      {
        hist[v = im[dy*xs+dxb]]++;
        if (v>vmed)
          nrig++;
        else
          if (v<vmed)
            nlef++;
      }
      // move median
      while (nrig>hist[vmed]+nlef)
      {
        nlef += hist[vmed];
        vmed++;
        nrig -= hist[vmed];
      }
      while (nlef>hist[vmed]+nrig)
      {
        nrig += hist[vmed];
        vmed--;
        nlef -= hist[vmed];
      }
      // save value
      line[i] = (unsigned char)vmed;
      // remove points from left
      dxb = MIA_bound(i-kew,0,xs-1);
      for (dy=-kew;dy<=kew;dy++)
      {
        hist[v = im[dy*xs+dxb]]--;
        if (v>vmed)
          nrig--;
        else
          if (v<vmed)
            nlef--;
      }
    }
    return ERR_OK;
}

// make median of all array for short
int IPL_HIST_mdn_pixS(
  const short* bufin, // input array
  int          lstr)  // size bufin
{                     // Returns: median
   int i,half,S;
   short asARM[512];

     MIA_memset(&asARM[0],0,sizeof(asARM));
     half = (lstr/2)|1;
     for (i=0;i<lstr;asARM[bufin[i++]&511]++);   // !!!MIA
     for (i=0,S=0;S<half;S+=asARM[i++]);
     return i-1;
}

// normalize histogram to a given value
MIA_RESULT_CODE IPL_HIST_Normalize(
  int* dst,         // OUT: normalized histogram
  const int* src,   // IN:  source histogram, if NULL, transform is in-place for dst
  int len,          // IN:  number of elements
  int maxabsval)    // IN:  maximum absolute value
{
  int maxval,minval,cIdx;

    // check arguments
    if (dst==NULL)
      return ERR_GEN_NULLPOINTER;
    if (src==NULL)
      src = dst;
    if (len<=0)
      return ERR_GEN_NO_DATA;
    if (maxabsval<=0)
      return ERR_GEN_INVALID_PARAMETER;
    // get max and min
    maxval = -0x7fffffff;
    minval =  0x7fffffff;
    for (cIdx=0;cIdx<len;cIdx++)
    {
      if (minval>src[cIdx])
        minval = src[cIdx];
      if (maxval<src[cIdx])
        maxval = src[cIdx];
    }
    // possibly negative has bigger abs value
    if (minval<0)
      minval = -minval;
    if (maxval<minval)
      maxval = minval;
    // check for zero
    if (maxval==0)
    {
      if (src!=dst)
        memset(dst,0,len*sizeof(dst[0]));
      return ERR_OK;
    }
    if (maxabsval!=maxval)
      // we need renormalization
      for (cIdx=0;cIdx<len;cIdx++)
        dst[cIdx] = (src[cIdx]*maxabsval+((src[cIdx]>=0)?maxval:-maxval)/2)/maxval;
    else
      // no renormalization
      if (src!=dst)
        // but still we need copying
        memcpy(dst,src,len*sizeof(dst[0]));
    return ERR_OK;
}

MIA_RESULT_CODE IPL_HIST_EqualizeHist(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys)
{
  int x,y,sum,i,val,hist[256];
  double scale;
  unsigned char lut[256+1];
  unsigned char* dptr;
  const unsigned char* sptr;

    // check arguments
    if ((dst==NULL)||(src==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((xs<=0)||(ys<=0))
      return ERR_GEN_NO_DATA;
    memset(hist, 0, sizeof(hist));
    for( y = 0; y < ys; y++ )
    {
      sptr = src+xs*y;
      for( x = 0; x < xs; x++ )
        hist[sptr[x]]++;
    }
    scale = 255.f/(xs*ys);
    sum = 0;
    for( i = 0; i < 256; i++ )
    {
      sum += hist[i];
      val = (int)(sum*scale+.5);
      lut[i] = (unsigned char)val;
    }
    lut[0] = 0;
    for( y = 0; y < ys; y++ )
    {
      sptr = src+xs*y;
      dptr = dst+xs*y;
      for( x = 0; x < xs; x++ )
        dptr[x] = lut[sptr[x]];
    }
    return ERR_OK;
}

// with strides
MIA_RESULT_CODE IPL_HIST_EqualizeHistExt(
  unsigned char* dst,
  int dststr,
  const unsigned char* src,
  int srcstr,
  int xs,
  int ys)
{
  int x,y,sum,i,val,hist[256];
  double scale;
  unsigned char lut[256+1];
  unsigned char* dptr;
  const unsigned char* sptr;

    // check arguments
    if ((dst==NULL)||(src==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((xs<=0)||(ys<=0))
      return ERR_GEN_NO_DATA;
    memset(hist, 0, sizeof(hist));
    for( y = 0; y < ys; y++ )
    {
      sptr = src+y*srcstr;
      for( x = 0; x < xs; x++ )
        hist[sptr[x]]++;
    }
    scale = 255.f/(xs*ys);
    sum = 0;
    for( i = 0; i < 256; i++ )
    {
      sum += hist[i];
      val = (int)(sum*scale+.5);
      lut[i] = (unsigned char)val;
    }
    lut[0] = 0;
    for( y = 0; y < ys; y++ )
    {
      sptr = src+srcstr*y;
      dptr = dst+dststr*y;
      for( x = 0; x < xs; x++ )
        dptr[x] = lut[sptr[x]];
    }
    return ERR_OK;
}
