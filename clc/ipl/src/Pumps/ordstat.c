/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  16 March 2006                                          */
/*      Revision: 1.0.00                                                 */
/*      Purpose:  order statistic filters - median, morphology, etc      */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/
#define __FILENUM__ 5 // __FILENUM__TAG5

#include <stdlib.h>
#include "../ownipl.h"

MIA_RESULT_CODE IPL_FILT_RFilter(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  int kehw,
  int kehh)
{
  int i,j,k,vmin,vmax,vmed,v,nlef,nrig,dir,j_fin;
  unsigned char hist[256];

    if ((dst==NULL)||(src==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((xs<=2*kehw)||(ys<=2*kehw)||(kehw<=0)||(kehh<=0))
      return ERR_GEN_INVALID_PARAMETER;
    // fill in the window histogram and calculate min,max,med
    vmin = 255;
    vmax = 0;
    memset(hist,0,sizeof(hist));
    for (i=0;i<2*kehh;i++)
      for (j=0;j<2*kehw;j++)
      {
        hist[v = src[i*xs+j]]++;
        if (vmin>v)
          vmin = v;
        if (vmax<v)
          vmax = v;
      }
    // median
    vmed = 0;
    nlef = 0;
    nrig = 4*kehh*kehw-hist[0];
    //
    dir = 1;
    j = kehw;
    // vertical movement
    for (i=kehh;i<ys-kehh;i++)
    {
      // vertical add
      for (k=-kehw;k<kehw;k++)
      {
        hist[v = src[(i+kehh)*xs+j+dir*k]]++;
        if (vmin>v)
          vmin = v;
        if (vmax<v)
          vmax = v;
        if (v>vmed)
          nrig++;
        else
          if (v<vmed)
            nlef++;
      }
      // horizontal movement
      j_fin=j+dir*(xs-2*kehw);
      while (j!=j_fin)
      {
        // horizontal add
        for (k=-kehh;k<=kehh;k++)
        {
          hist[v = src[(i+k)*xs+j+dir*kehw]]++;
          if (vmin>v)
            vmin = v;
          if (vmax<v)
            vmax = v;
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
        // calc destination
        dst[i*xs+j] = (unsigned char)((vmed>(vmax+vmin)/2)?vmax:((vmed<(vmax+vmin)/2)?vmin:(vmax+vmin)/2));
//        dst[i*xs+j] = (unsigned char)((vmed>127)?vmax:((vmed<127)?vmin:(vmax+vmin)/2));
//        dst[i*xs+j] = (unsigned char)((vmed>(vmax+vmin)/2)?255:0);
//        dst[i*xs+j] = (unsigned char)(vmax-vmin);
        // horizontal sub
        for (k=-kehh;k<=kehh;k++)
        {
          if (!(--hist[v = src[(i+k)*xs+j-dir*kehw]]))
          {
            if (vmin==v)
              while (!hist[vmin])
                vmin++;
            else
              if (vmax==v)
                while (!hist[vmax])
                  vmax--;
          }
          if (v>vmed)
            nrig--;
          else
            if (v<vmed)
              nlef--;
        }
        // move pointer
        j += dir;
      }
      // vertical sub
      for (k=-kehw;k<kehw;k++)
      {
        if (!(--hist[v = src[(i-kehh)*xs+j+dir*k]]))
        {
          if (vmin==v)
            while (!hist[vmin])
              vmin++;
          else
            if (vmax==v)
              while (!hist[vmax])
                vmax--;
        }
        if (v>vmed)
          nrig--;
        else
          if (v<vmed)
            nlef--;
      }
      // refine horizontal ptr
      j -= dir;
      // change direction
      dir = -dir;
    }
    // fill borders
    ownIPL_FillBordersXB(dst,xs,ys,kehw,kehh);
    return ERR_OK;
}

void IPL_FILT_Median(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  int kehw,
  int kehh)
{
  int i,j,k,vmed,v,nlef,nrig,dir,j_fin;
  unsigned char hist[256];

    // fill in the window histogram and calculate min,max,med
    memset(hist,0,sizeof(hist));
    for (i=0;i<2*kehh;i++)
      for (j=0;j<2*kehw;j++)
        hist[v = src[i*xs+j]]++;
    // median
    vmed = 0;
    nlef = 0;
    nrig = 4*kehh*kehw-hist[0];
    //
    dir = 1;
    j = kehw;
    // vertical movement
    for (i=kehh;i<ys-kehh;i++)
    {
      // vertical add
      for (k=-kehw;k<kehw;k++)
      {
        hist[v = src[(i+kehh)*xs+j+dir*k]]++;
        if (v>vmed)
          nrig++;
        else
          if (v<vmed)
            nlef++;
      }
      // horizontal movement
      j_fin=j+dir*(xs-2*kehw);
      while (j!=j_fin)
      {
        // horizontal add
        for (k=-kehh;k<=kehh;k++)
        {
          hist[v = src[(i+k)*xs+j+dir*kehw]]++;
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
        // calc destination
        dst[i*xs+j] = (unsigned char)vmed;
        // horizontal sub
        for (k=-kehh;k<=kehh;k++)
        {
          hist[v = src[(i+k)*xs+j-dir*kehw]]--;
          if (v>vmed)
            nrig--;
          else
            if (v<vmed)
              nlef--;
        }
        // move pointer
        j += dir;
      }
      // vertical sub
      for (k=-kehw;k<kehw;k++)
      {
        hist[v = src[(i-kehh)*xs+j+dir*k]]--;
        if (v>vmed)
          nrig--;
        else
          if (v<vmed)
            nlef--;
      }
      // refine horizontal ptr
      j -= dir;
      // change direction
      dir = -dir;
    }
    // fill borders
    ownIPL_FillBordersXB(dst,xs,ys,kehw,kehh);
}

void IPL_FILT_Dilate(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  int kehw,
  int kehh)
{
  int i,j,k,vmax,v,dir,j_fin,hist[256];

    // fill in the window histogram and calculate min,max,med
    memset(hist,0,sizeof(hist));
    vmax = 0;
    for (i=0;i<2*kehh;i++)
      for (j=0;j<2*kehw;j++)
      {
        hist[v = src[i*xs+j]]++;
        if (vmax<v)
          vmax = v;
      }
    //
    dir = 1;
    j = kehw;
    // vertical movement
    for (i=kehh;i<ys-kehh;i++)
    {
      // vertical add
      for (k=-kehw;k<kehw;k++)
      {
        hist[v = src[(i+kehh)*xs+j+dir*k]]++;
        if (vmax<v)
          vmax = v;
      }
      // horizontal movement
      j_fin=j+dir*(xs-2*kehw);
      while (j!=j_fin)
      {
        // horizontal add
        for (k=-kehh;k<=kehh;k++)
        {
          hist[v = src[(i+k)*xs+j+dir*kehw]]++;
          if (vmax<v)
            vmax = v;
        }
        // set destination
        dst[i*xs+j] = (unsigned char)vmax;
        // horizontal sub
        for (k=-kehh;k<=kehh;k++)
          hist[src[(i+k)*xs+j-dir*kehw]]--;
        // move max
        while (hist[vmax]==0)
          vmax--;
        // move pointer
        j += dir;
      }
      // vertical sub
      for (k=-kehw;k<kehw;k++)
        hist[src[(i-kehh)*xs+j+dir*k]]--;
      // move max
      while (hist[vmax]==0)
        vmax--;
      // refine horizontal ptr
      j -= dir;
      // change direction
      dir = -dir;
    }
    // fill borders
    ownIPL_FillBordersXB(dst,xs,ys,kehw,kehh);
}

void IPL_FILT_Erode(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys,
  int kehw,
  int kehh)
{
  int i,j,k,vmin,v,dir,j_fin,hist[256];

    // fill in the window histogram and calculate min,max,med
    memset(hist,0,sizeof(hist));
    vmin = 255;
    for (i=0;i<2*kehh;i++)
      for (j=0;j<2*kehw;j++)
      {
        hist[v = src[i*xs+j]]++;
        if (vmin>v)
          vmin = v;
      }
    //
    dir = 1;
    j = kehw;
    // vertical movement
    for (i=kehh;i<ys-kehh;i++)
    {
      // vertical add
      for (k=-kehw;k<kehw;k++)
      {
        hist[v = src[(i+kehh)*xs+j+dir*k]]++;
        if (vmin>v)
          vmin = v;
      }
      // horizontal movement
      j_fin=j+dir*(xs-2*kehw);
      while (j!=j_fin)
      {
        // horizontal add
        for (k=-kehh;k<=kehh;k++)
        {
          hist[v = src[(i+k)*xs+j+dir*kehw]]++;
          if (vmin>v)
            vmin = v;
        }
        // set destination
        dst[i*xs+j] = (unsigned char)vmin;
        // horizontal sub
        for (k=-kehh;k<=kehh;k++)
          hist[src[(i+k)*xs+j-dir*kehw]]--;
        // move max
        while (hist[vmin]==0)
          vmin++;
        // move pointer
        j += dir;
      }
      // vertical sub
      for (k=-kehw;k<kehw;k++)
        hist[src[(i-kehh)*xs+j+dir*k]]--;
      // move max
      while (hist[vmin]==0)
        vmin++;
      // refine horizontal ptr
      j -= dir;
      // change direction
      dir = -dir;
    }
    // fill borders
    ownIPL_FillBordersXB(dst,xs,ys,kehw,kehh);
}

// 1D erosion for int type
void IPL_FILT_Erode1D_Cycled_Int(
  int* dst,         // OUT: destination array
  const int* src,   // IN:  source array
  int len,          // IN:  size of array
  int kehw)         // IN:  half aperture
{
  int i,j,minval;

    for (i=0;i<len;i++)
    {
      minval = src[(len+i-kehw)%len];
      for (j=-kehw+1;j<=kehw;j++)
        if (minval>src[(len+i+j)%len])
          minval = src[(len+i+j)%len];
      dst[i] = minval;
    }
}

// 1D dilation for int type
void IPL_FILT_Dilate1D_Cycled_Int(
  int* dst,         // OUT: destination array
  const int* src,   // IN:  source array
  int len,          // IN:  size of array
  int kehw)         // IN:  half aperture
{
  int i,j,maxval;

    for (i=0;i<len;i++)
    {
      maxval = src[(len+i-kehw)%len];
      for (j=-kehw+1;j<=kehw;j++)
        if (maxval<src[(len+i+j)%len])
          maxval = src[(len+i+j)%len];
      dst[i] = maxval;
    }
}

void IPL_FILT_Erode3x3(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys)
{
  int i,j,m,q;

    for (i=1;i<ys-1;i++)
    {
      for (j=1;j<xs-1;j++)
      {
        m =        src[(i-1)*xs+(j-1)];
        if (m>(q = src[(i-1)*xs+(j  )]))
          m = q;
        if (m>(q = src[(i-1)*xs+(j+1)]))
          m = q;
        if (m>(q = src[(i  )*xs+(j-1)]))
          m = q;
        if (m>(q = src[(i  )*xs+(j  )]))
          m = q;
        if (m>(q = src[(i  )*xs+(j+1)]))
          m = q;
        if (m>(q = src[(i+1)*xs+(j-1)]))
          m = q;
        if (m>(q = src[(i+1)*xs+(j  )]))
          m = q;
        if (m>(q = src[(i+1)*xs+(j+1)]))
          m = q;
        dst[i*xs+j] = (unsigned char)m;
      }
    }
    ownIPL_FillBorders1B(dst,xs,ys);
}

void IPL_FILT_Dilate3x3(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys)
{
  int i,j,m,q;

    for (i=1;i<ys-1;i++)
    {
      for (j=1;j<xs-1;j++)
      {
        m =        src[(i-1)*xs+(j-1)];
        if (m<(q = src[(i-1)*xs+(j  )]))
          m = q;
        if (m<(q = src[(i-1)*xs+(j+1)]))
          m = q;
        if (m<(q = src[(i  )*xs+(j-1)]))
          m = q;
        if (m<(q = src[(i  )*xs+(j  )]))
          m = q;
        if (m<(q = src[(i  )*xs+(j+1)]))
          m = q;
        if (m<(q = src[(i+1)*xs+(j-1)]))
          m = q;
        if (m<(q = src[(i+1)*xs+(j  )]))
          m = q;
        if (m<(q = src[(i+1)*xs+(j+1)]))
          m = q;
        dst[i*xs+j] = (unsigned char)m;
      }
    }
    ownIPL_FillBorders1B(dst,xs,ys);
}

void IPL_FILT_Erode3x3Cross(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys)
{
  int i,j,m,q;

    for (i=1;i<ys-1;i++)
    {
      for (j=1;j<xs-1;j++)
      {
        m =        src[(i-1)*xs+(j  )];
        if (m>(q = src[(i  )*xs+(j-1)]))
          m = q;
        if (m>(q = src[(i  )*xs+(j  )]))
          m = q;
        if (m>(q = src[(i  )*xs+(j+1)]))
          m = q;
        if (m>(q = src[(i+1)*xs+(j  )]))
          m = q;
        dst[i*xs+j] = (unsigned char)m;
      }
    }
    ownIPL_FillBorders1B(dst,xs,ys);
}

void IPL_FILT_Dilate3x3Cross(
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys)
{
  int i,j,m,q;

    for (i=1;i<ys-1;i++)
    {
      for (j=1;j<xs-1;j++)
      {
        m =  src[(i-1)*xs+(j  )];

        if (m<(q = src[(i  )*xs+(j-1)]))
          m = q;
        if (m<(q = src[(i  )*xs+(j  )]))
          m = q;
        if (m<(q = src[(i  )*xs+(j+1)]))
          m = q;
        if (m<(q = src[(i+1)*xs+(j  )]))
          m = q;
        dst[i*xs+j] = (unsigned char)m;
      }
    }
    ownIPL_FillBorders1B(dst,xs,ys);
}

void IPL_FILT_Erode5x1( //  along x
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys)
{
  int i,j,m,q;

    for (i=0;i<ys;i++)
    {
      for (j=2;j<xs-2;j++)
      {
        m =        src[i*xs+(j-2)];
        if (m>(q = src[i*xs+(j-1)]))
          m = q;
        if (m>(q = src[i*xs+(j  )]))
          m = q;
        if (m>(q = src[i*xs+(j+1)]))
          m = q;
        if (m>(q = src[i*xs+(j+2)]))
          m = q;
        dst[i*xs+j] = (unsigned char)m;
      }
    }
    for (i=0;i<ys;i++)
    {
      dst[i*xs+0] = dst[i*xs+1] = dst[i*xs+2];
      dst[i*xs+xs-1] = dst[i*xs+xs-2] = dst[i*xs+xs-3];
    }
}

void IPL_FILT_Dilate5x1( //  along x
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys)
{
  int i,j,m,q;

    for (i=0;i<ys;i++)
    {
      for (j=2;j<xs-2;j++)
      {
        m =        src[i*xs+(j-2)];
        if (m<(q = src[i*xs+(j-1)]))
          m = q;
        if (m<(q = src[i*xs+(j  )]))
          m = q;
        if (m<(q = src[i*xs+(j+1)]))
          m = q;
        if (m<(q = src[i*xs+(j+2)]))
          m = q;
        dst[i*xs+j] = (unsigned char)m;
      }
    }
    for (i=0;i<ys;i++)
    {
      dst[i*xs+0] = dst[i*xs+1] = dst[i*xs+2];
      dst[i*xs+xs-1] = dst[i*xs+xs-2] = dst[i*xs+xs-3];
    }
}

void IPL_FILT_Erode1x5( //  along y
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys)
{
  int i,j,m,q;

    for (i=2;i<ys-2;i++)
    {
      for (j=0;j<xs;j++)
      {
        m =        src[(i-2)*xs+j];
        if (m>(q = src[(i-1)*xs+j]))
          m = q;
        if (m>(q = src[(i  )*xs+j]))
          m = q;
        if (m>(q = src[(i+1)*xs+j]))
          m = q;
        if (m>(q = src[(i+2)*xs+j]))
          m = q;
        dst[i*xs+j] = (unsigned char)m;
      }
    }
    memcpy(dst,dst+xs*2,xs);
    memcpy(dst+xs,dst+xs*2,xs);
    memcpy(dst+xs*(ys-2),dst+xs*(ys-3),xs);
    memcpy(dst+xs*(ys-1),dst+xs*(ys-3),xs);
}

void IPL_FILT_Dilate1x5( //  along y
  unsigned char* dst,
  const unsigned char* src,
  int xs,
  int ys)
{
  int i,j,m,q;

    for (i=2;i<ys-2;i++)
    {
      for (j=0;j<xs;j++)
      {
        m =        src[(i-2)*xs+j];
        if (m<(q = src[(i-1)*xs+j]))
          m = q;
        if (m<(q = src[(i  )*xs+j]))
          m = q;
        if (m<(q = src[(i+1)*xs+j]))
          m = q;
        if (m<(q = src[(i+2)*xs+j]))
          m = q;
        dst[i*xs+j] = (unsigned char)m;
      }
    }
    memcpy(dst,dst+xs*2,xs);
    memcpy(dst+xs,dst+xs*2,xs);
    memcpy(dst+xs*(ys-2),dst+xs*(ys-3),xs);
    memcpy(dst+xs*(ys-1),dst+xs*(ys-3),xs);
}

#define PIX_SORT(x,y) v=(b##x<b##y)?b##x:b##y; b##y+=b##x-v; b##x=v;

MIA_RESULT_CODE IPL_FILT_Median3x3(
  unsigned char* dst,           // destination
  const unsigned char* src,     // source
  int xs,               // size of image in elements
  int ys)
{
  int i,j;
  unsigned int b0,b1,b2,b3,b4,b5,b6,b7,b8,v;

    src += xs+1;
    for (i=1;i<ys-1;i++)
    {
      for (j=1;j<xs-1;j++)
      {
        b0 = src[-xs-1];
        b1 = src[-xs  ];
        b2 = src[-xs+1];
        b3 = src[   -1];
        b4 = src[    0];
        b5 = src[    1];
        b6 = src[ xs-1];
        b7 = src[ xs  ];
        b8 = src[ xs+1];
        PIX_SORT(1,2);
        PIX_SORT(4,5);
        PIX_SORT(7,8);
        PIX_SORT(0,1);
        PIX_SORT(3,4);
        PIX_SORT(6,7);
        PIX_SORT(1,2);
        PIX_SORT(4,5);
        PIX_SORT(7,8);
        PIX_SORT(0,3);
        PIX_SORT(5,8);
        PIX_SORT(4,7);
        PIX_SORT(3,6);
        PIX_SORT(1,4);
        PIX_SORT(2,5);
        PIX_SORT(4,7);
        PIX_SORT(4,2);
        PIX_SORT(6,4);
        PIX_SORT(4,2);
        dst[i*xs+j] = (unsigned char)b4;
        src++;
      }
      src += 2;
    }
    // fill left and right margins
    for (i=1;i<ys-1;i++)
    {
      dst[i*xs] = dst[i*xs+1];
      dst[i*xs+xs-1] = dst[i*xs+xs-2];
    }
    MIA_memcpy(dst,dst+xs,xs);
    MIA_memcpy(dst+(ys-1)*xs,dst+(ys-2)*xs,xs);
    return ERR_OK;
}

#define SIGN(I) ((I)<0?0:(I)==0?1:2)

void IPL_FILT_mdn2(
  const unsigned char *bufin,
  unsigned char *bufout,
  int lstr,
  int apt_size)
{
  int i,j,k;
  int half1_apt,half2_apt;
  int MED,S,N_SLEFT;
  int v_old,v_new;
  unsigned char ARM[256];

    half1_apt=apt_size/2;
    half2_apt=apt_size-half1_apt-1;
    MIA_memset(&ARM[0],0,sizeof(ARM));
    for (i=0;i<half1_apt;i++) 
      ARM[bufin[lstr-1-i]]++;
    for (i=half1_apt,j=0;i<apt_size;i++,j++) 
      ARM[bufin[j]]++;
    for (i=0,S=0;S<half1_apt+1;) S+=ARM[i++];
    MED=i-1;
    N_SLEFT=S;
    bufout[0]=(unsigned char)MED;
    for (i=1;i<lstr;i++)
    {
      k=i-half1_apt-1;
      v_old=(k<0)?bufin[lstr-half1_apt-1+i]:bufin[k];
      k=i+half2_apt;
      v_new=(k>=lstr)?bufin[lstr-i-1]:bufin[k];
      ARM[v_old]--;
      ARM[v_new]++;
      if (SIGN(v_old-MED)==SIGN(v_new-MED))
      {
        bufout[i]=(unsigned char)MED;
        continue;
      }
      else if (v_new==MED)
      {
        bufout[i]=(unsigned char)MED;
        if (v_old>MED)
          N_SLEFT++;
        continue;
      }
      if (v_new>MED)
      {
        N_SLEFT--;
        if (N_SLEFT<(half1_apt+1))
        {
          for (j=MED+1;ARM[j++]==0;);
          j--;
          MED=j;
          N_SLEFT+=ARM[j];
        }
        bufout[i]=(unsigned char)MED;
      }
      else
      {
        if (v_old!=MED)
          N_SLEFT++;
        if ((N_SLEFT-ARM[MED])>(half1_apt+1))
        {
          N_SLEFT-=ARM[MED];
          for (j=MED-1;ARM[j--]==0;);
          MED=j+1;
        }
        bufout[i]=(unsigned char)MED;
      }
    }
}

void IPL_FILT_mdnS256(
  short *bufin,
  short *bufout,
  int lstr,
  int apt_size)
{
  int i,j,k;
  int half1_apt,half2_apt;
  int MED,S,N_SLEFT;
  int v_old,v_new;
  unsigned char ARM[256];

    half1_apt=apt_size/2;
    half2_apt=apt_size-half1_apt-1;
    MIA_memset(&ARM[0],0,sizeof(ARM));
    for (i=0;i<half1_apt;i++) 
      ARM[bufin[lstr-1-i]]++;
    for (i=half1_apt,j=0;i<apt_size;i++,j++) 
      ARM[bufin[j]]++;
    for (i=0,S=0;S<half1_apt+1;) S+=ARM[i++];
    MED=i-1;
    N_SLEFT=S;
    bufout[0]=(short)MED;
    for (i=1;i<lstr;i++)
    {
      k=i-half1_apt-1;
      v_old=(k<0)?bufin[lstr-half1_apt-1+i]:bufin[k];
      k=i+half2_apt;
      v_new=(k>=lstr)?bufin[lstr-i-1]:bufin[k];
      ARM[v_old]--;
      ARM[v_new]++;
      if (SIGN(v_old-MED)==SIGN(v_new-MED))
      {
        bufout[i]=(short)MED;
        continue;
      }
      else if (v_new==MED)
      {
        bufout[i]=(short)MED;
        if (v_old>MED)
          N_SLEFT++;
        continue;
      }
      if (v_new>MED)
      {
        N_SLEFT--;
        if (N_SLEFT<(half1_apt+1))
        {
          for (j=MED+1;ARM[j++]==0;);
          j--;
          MED=j;
          N_SLEFT+=ARM[j];
        }
        bufout[i]=(short)MED;
      }
      else
      {
        if (v_old!=MED)
          N_SLEFT++;
        if ((N_SLEFT-ARM[MED])>(half1_apt+1))
        {
          N_SLEFT-=ARM[MED];
          for (j=MED-1;ARM[j--]==0;);
          MED=j+1;
        }
        bufout[i]=(short)MED;
      }
    }
}

staticFunc int ownIPL_sortint_linemed(
  const void* e1,
  const void* e2)
{
    return ((int*)e1)[0]-((int*)e2)[0];
}

staticFunc int ownIPL_sortdouble_linemed(
  const void* e1,
  const void* e2)
{
    if (((double*)e1)[0]>((double*)e2)[0])
      return 1;
    if (((double*)e1)[0]<((double*)e2)[0])
      return -1;
    return 0;
}

// straight-thru median for sequence of int's
MIA_RESULT_CODE IPL_FILT_LineMedian_Int(
  int *dst,         // OUT: output array
  const int* src,   // IN:  input array
  int len,          // IN:  length of input
  int hw)           // IN:  half-width of median window
{
  int arra[100],cIdx,cWin;

    // check arguments
    if ((src==NULL)||(dst==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((len<=0)||(hw<0))
      return ERR_GEN_NO_DATA;
    if (2*hw+1>sizeof(arra)/sizeof(arra[0]))
      return ERR_GEN_INSUFFICIENT_BUFFER;
    // process
    for (cIdx=0;cIdx<len;cIdx++)
    {
      for (cWin=-hw;cWin<=hw;cWin++)
        if (cIdx+cWin<0)
          arra[cWin+hw] = src[0];
        else
          if (cIdx+cWin>=len)
            arra[cWin+hw] = src[len-1];
          else
            arra[cWin+hw] = src[cIdx+cWin];
      qsort(arra,hw*2+1,sizeof(arra[0]),ownIPL_sortint_linemed);
      dst[cIdx] = arra[hw];
    }
    return ERR_OK;
}

// straight-thru median for sequence of double's
MIA_RESULT_CODE IPL_FILT_LineMedian_Double(
  double *dst,        // OUT: output array
  const double* src,  // IN:  input array
  int len,            // IN:  length of input
  int hw)             // IN:  half-width of median window
{
  double arra[100];
  int cIdx,cWin;

    // check arguments
    if ((src==NULL)||(dst==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((len<=0)||(hw<0))
      return ERR_GEN_NO_DATA;
    if (2*hw+1>sizeof(arra)/sizeof(arra[0]))
      return ERR_GEN_INSUFFICIENT_BUFFER;
    // process
    for (cIdx=0;cIdx<len;cIdx++)
    {
      for (cWin=-hw;cWin<=hw;cWin++)
        if (cIdx+cWin<0)
          arra[cWin+hw] = src[0];
        else
          if (cIdx+cWin>=len)
            arra[cWin+hw] = src[len-1];
          else
            arra[cWin+hw] = src[cIdx+cWin];
      qsort(arra,hw*2+1,sizeof(arra[0]),ownIPL_sortdouble_linemed);
      dst[cIdx] = arra[hw];
    }
    return ERR_OK;
}

// Dilated image should be binary, i.e. contain only two values: 
// zero and non-zero (and ALL non-zero values should be equal)
// 255 is used in dst as non-zero value. 
MIA_RESULT_CODE IPL_FILT_Dilate_Binary(
  unsigned char* dst,         // OUT: resulting image
  const unsigned char* src,   // IN:  original image
  int xs,                     // IN:  image size
  int ys,                     //
  int kehw,                   // IN:  kernel half-size
  int kehh,                   //
  void* pvBuf,                // IN:  temporary buffer
  int* pnLen)                 // IN/OUT: bytes allocated/used
{
//  const unsigned char* psrc;
  int *sums,sz,i,j,Sum;

    // check arguments
    if (pvBuf==NULL)
      return ERR_GEN_NULLPOINTER;
    if ((kehw<=0)||(kehh<=0))
      return ERR_GEN_NO_DATA;
    if ((2*kehw>=xs)||(2*kehh>=ys))
      return ERR_GEN_INVALID_PARAMETER;
    // calc size
    sz = sizeof(sums[0])*xs;
    if (pvBuf==NULL)
    { // this is mem query
      *pnLen = sz;
      return ERR_OK;
    }
    // more check arguments
    if ((src==NULL)||(dst==NULL))
      return ERR_GEN_NULLPOINTER;
    if (*pnLen<sz)
    { // this is mem query
      *pnLen = sz;
      return ERR_GEN_INSUFFICIENT_BUFFER;
    }
    *pnLen = sz;
    // allocate
    sums = (int*)pvBuf;
    // init sums
/*    for (j=xs;j;j--)
    {
      *sums = 0;
      for (i=2*kehh;i;i--)
      {
        *sums += *src;
        src += xs;
      }
      sums++;
      src -= 2*kehh*xs-1;
    }
    sums -= xs;
    src -= xs;*/
    for (j=0;j<xs;j++)
      sums[j] = 0;
    for (i=0;i<2*kehh*xs;i+=xs)
      for (j=0;j<xs;j++)
        sums[j] += src[i+j];
    // process (central part)
    dst += kehh*xs;
    for (i=ys-2*kehh;i;i--)
    {
      // add bottom row to sums
      src += (2*kehh)*xs;
      for (j=0;j<xs;j++)
        sums[j] += src[j];
      src -= (2*kehh)*xs;
      // init Sum
      Sum = 0;
      for (j=0;j<2*kehw;j++)
        Sum += sums[j];
      for (j=kehw;j<xs-kehw;j++)
      {
        // add right sum
        Sum += sums[j+kehw];
        // fill in dst
        dst[j] = (Sum)?255:0;
        // sub left sum
        Sum -= sums[j-kehw];
      }
      // sub top row from sums
      for (j=0;j<xs;j++)
        sums[j] -= src[j];
      // shift pointers
      src += xs;
      dst += xs;
    }
    dst -= (ys-kehh)*xs;
    // fill borders
    ownIPL_FillBordersXB(dst,xs,ys,kehw,kehh);
    return ERR_OK;
}

void IPL_Morph_DilateLineBinaryCycled(
  int* dst,
  const int* src,
  int len,
  int wid)
{
  int i,j;

    for(i=0;i<len;i++)
    {
      for (j=-wid;j<=wid;j++)
        if (src[(i+j+len)%len])
          break;
      dst[i] = (j>wid)?0:1;
    }
}

void IPL_Morph_ErodeLineBinaryCycled(
  int* dst,
  const int* src,
  int len,
  int wid)
{
  int i,j;

    for(i=0;i<len;i++)
    {
      for (j=-wid;j<=wid;j++)
        if (!src[(i+j+len)%len])
          break;
      dst[i] = (j>wid)?1:0;
    }
}
