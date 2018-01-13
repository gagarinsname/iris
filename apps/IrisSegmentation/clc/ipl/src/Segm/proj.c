/*----------------------------------------------------------------------------*/
/*                                                                            */
/*      Created:  4 August 2010                                               */
/*      Revision: 1.0.00                                                      */
/*      Purpose:  projections                                                 */
/*      Authors:                                                              */
/*        Ivan Matveev                                                        */
/*                                                                            */
/*----------------------------------------------------------------------------*/
#define __FILENUM__ 9 // __FILENUM__TAG9

#include <stdio.h>
#include <malloc.h>
#include <memory.h>
#include "bpl.h"
#include "ipl.h"
//#include "dbg.h"

#define BLUR_HIST_WND 4

//#define BOTH_GRAD_DIR
#undef BOTH_GRAD_DIR

// with blurring
MIA_RESULT_CODE IPL_PROJ_FindHoughDonatorProjection7(
  int* pnProjRN_b,          // OUT: circular projection - right side nums
  int* pnProjLN_b,          // OUT: circular projection - left side nums
  const unsigned char* im,  // IN:  image
  int xs,                   // IN:  image size
  int ys,                   // 
  int xc,                   // IN:  ring center
  int yc,                   //
  int r,                    // IN:  inner radius
  int* pR,                  // IN/OUT: outer radius proposed/actually used
  int nMaxRadOfDence,       // IN:  maximum radius for dense processing. 0 - dense always
  void* buf,                // IN:  external buffer
  int* buflen)              // IN/OUT: allocated/used bytes
{
  int i,j,gx,gy,y,G,L,S,*numpoR,*numpoL,currad,R;
  int j_inn,j_out,Y,D,RR_ii,rr_ii;
  int *pnProjRN,*pnProjLN,*numpoR_b,*numpoL_b;
  short *smoothedim,*psim;
  double* tempim;

    // check arguments
    if ((buflen==NULL)||(pR==NULL))
      return ERR_GEN_NULLPOINTER;
    // verify radii and calculate max outer radius
    R = *pR;
    if (R<10)
      return ERR_GEN_SIZE_NOT_MATCH;
    if ((r+5>R)||(r<5))
      return ERR_GEN_SIZE_NOT_MATCH;
    if ((xs<2*r)||(ys<2*r))
      return ERR_GEN_NO_DATA;
    if ((xc-r/2<0)||(xc+r/2>=xs)||(yc-r/2<0)||(yc+r/2>=ys))
      return ERR_GEN_INVALID_PARAMETER;
    // verify outer radius - bound by maximum span from center to border
    i = xs-1-xc;
    if (i<xc)
      i = xc;
    if (R>i)
      R = i;
    // verify outer radius - bound by image size
    if (2*R>xs)
      R = xs/2;
    // return corrected outer radius
    *pR = R;
    // calc size
    i = (2+4+2)*(R+1)*sizeof(numpoR[0])+
        xs*ys*sizeof(smoothedim[0])+xs*ys*sizeof(tempim[0]);
    if (buf==NULL)
    {
      *buflen = i;
      return ERR_OK;
    }
    if (*buflen<i)
    {
      *buflen = i;
      return ERR_GEN_NOMEMORY;
    }
    // more check arguments
    if ((pnProjRN_b==NULL)||(pnProjLN_b==NULL)||(im==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((((ptr_t)im)&3)||(((ptr_t)buf)&3)||(xs&3)||
        (((ptr_t)pnProjRN_b)&3)||(((ptr_t)pnProjLN_b)&3))
      return ERR_GEN_BAD_ALIGNMENT;
    // allocate
    tempim = (double*)buf;
    numpoR   = (int*)(tempim+xs*ys);
    numpoL   = numpoR+(R+1);
    pnProjRN = numpoL+(R+1);
    pnProjLN = pnProjRN+(R+1);
    numpoR_b = pnProjLN+(R+1);
    numpoL_b = numpoR_b+(R+1);
    smoothedim = (short*)(numpoL_b+(R+1));
    // clear junk
    MIA_memset(pnProjRN,0,sizeof(pnProjRN[0])*(1+R));
    MIA_memset(pnProjLN,0,sizeof(pnProjLN[0])*(1+R));
    MIA_memset(numpoR,0,sizeof(numpoR[0])*(1+R));
    MIA_memset(numpoL,0,sizeof(numpoL[0])*(1+R));
    IPL_FILT_GaussianSmooth_uint8(
      smoothedim,
      tempim,
      im,
      xs,
      ys,
      2.);
    // along y
    if (!nMaxRadOfDence)
      D = 0x10000;
    else
    {
      D = (r*0x10000)/nMaxRadOfDence;
      if (D<0x10000)
        D = 0x10000;
    }
    for (Y=-R*0x10000;Y<=R*0x10000;Y+=D)
    {
      i = (Y+0x8000)/0x10000;
      y = yc+i;
      if ((y<1)||(y>=ys-1))
        continue;
      // left side
      j_out = RR_ii = MIA_INT_sqrt(R*R-i*i);
      if (xc-j_out<1)
        j_out = xc-1;
      j_inn = (i>0)?i:(-i);
//j_inn = 2*((i>0)?i:(-i));
      if (j_inn<r)
      {
        rr_ii = MIA_INT_sqrt(r*r-i*i);
        if (j_inn<rr_ii)
          j_inn = rr_ii;
      }
      L = i*i+j_out*j_out;
      currad = MIA_INT_sqrt(L);
      psim = smoothedim+y*xs+xc-j_out;
      for (j=-j_out;j<=-j_inn;j++)
      {
        // count number of points in circle
        numpoL[currad]++;
        // straight gradient
//        gx = (int)(1.2*(psim[ 1]-psim[-1]));
        gx = 6*(psim[ 1]-psim[-1])/5;
        if (gx<-2)
        {
          gy =psim[ xs]-psim[-xs];
          G = gx*gx+gy*gy;
          if (G>2*36*2)
          {
            // inner product
            S = gx*j+gy*i;
            //
#ifdef BOTH_GRAD_DIR
            if (((int40)(S*3))*S>((int40)L)*G*2)
              pnProjLN[currad]++;
#else //BOTH_GRAD_DIR
//            if ((S>0)&&(((int40)(S*3))*S>((int40)L)*G*2))
if ((S>0)&&(((int40)(S*4))*S>((int40)L)*G*3))
//if ((S>0)&&(((int40)(S*5))*S>((int40)L)*G*4))
              pnProjLN[currad]++;
#endif //BOTH_GRAD_DIR
          }
        }
        psim++;
        L -= (-j*2-1);
        if (L<currad*currad)
          currad--;
      }
      // right side
      j_out = RR_ii;
      if (xc+j_out>xs-2)
        j_out = xs-2-xc;
      L = i*i+j_out*j_out;
      currad = MIA_INT_sqrt(L);
      psim = smoothedim+y*xs+xc+j_out;
      for (j=j_out;j>=j_inn;j--)
      {
        // count number of points in circle
        numpoR[currad]++;
        // straight gradient
//        gx = (int)(1.2*(psim[ 1]-psim[-1]));
        gx = 6*(psim[ 1]-psim[-1])/5;
        if (gx>2)
        {
          gy =psim[ xs]-psim[-xs];
          G = gx*gx+gy*gy;
          if (G>2*36*2)
          {
            // inner product
            S = gx*j+gy*i;
            //
#ifdef BOTH_GRAD_DIR
            if (((int40)(S*3))*S>((int40)L)*G*2)
              pnProjRN[currad]++;
#else //BOTH_GRAD_DIR
            //if ((S>0)&&(((int40)(S*3))*S>((int40)L)*G*2))
if ((S>0)&&(((int40)(S*4))*S>((int40)L)*G*3))
//if ((S>0)&&(((int40)(S*5))*S>((int40)L)*G*4))
              pnProjRN[currad]++;
#endif //BOTH_GRAD_DIR
          }
        }
        psim--;
        L -= (j*2-1);
        if (L<currad*currad)
          currad--;
      }
    }
    // normalize
    for (i=0;i<r;i++)
    {
      numpoR[i] = numpoR[r];
      numpoL[i] = numpoL[r];
      pnProjRN[i] = pnProjRN[r];
      pnProjLN[i] = pnProjLN[r];
    }
    for (i=0;i<=R;i++)
    {
      pnProjRN[i] = (numpoR[i])?((pnProjRN[i]*1024)/numpoR[i]):0;
      pnProjLN[i] = (numpoL[i])?((pnProjLN[i]*1024)/numpoL[i]):0;
    }
    // process histogram: blur
    IPL_HIST_Blur(numpoR_b,numpoR,R+1,BLUR_HIST_WND);
    IPL_HIST_Blur(numpoL_b,numpoL,R+1,BLUR_HIST_WND);
    IPL_HIST_Blur(pnProjRN_b,pnProjRN,R+1,BLUR_HIST_WND);
    IPL_HIST_Blur(pnProjLN_b,pnProjLN,R+1,BLUR_HIST_WND);
    // correct for small presence
    for (i=1;i<=R;i++)
    {
      L = 10000*numpoR_b[i]/(157*(2*BLUR_HIST_WND+1)*i);
      if (L<20)
        pnProjRN_b[i] = 0;
      L = 10000*numpoL_b[i]/(157*(2*BLUR_HIST_WND+1)*i);
      if (L<20)
        pnProjLN_b[i] = 0;
    }
    pnProjRN_b[0] = pnProjRN_b[1];
    pnProjLN_b[0] = pnProjLN_b[1];
    return ERR_OK;
}

// same as '5', but with thresholds from '3'
MIA_RESULT_CODE IPL_PROJ_FindHoughDonatorProjection6(
  int* pnProjRN_b,          // OUT: circular projection - right side nums
  int* pnProjLN_b,          // OUT: circular projection - left side nums
  const unsigned char* im,  // IN:  image
  int xs,                   // IN:  image size
  int ys,                   // 
  int xc,                   // IN:  ring center
  int yc,                   //
  int r,                    // IN:  inner radius
  int* pR,                  // IN/OUT: outer radius proposed/actually used
  int nMaxRadOfDence,       // IN:  maximum radius for dense processing. 0 - dense always
  void* buf,                // IN:  external buffer
  int* buflen)              // IN/OUT: allocated/used bytes
{
  int i,j,gx,gy,y,G,L,S,*numpoR,*numpoL,currad,R;
  int j_inn,j_out,Y,D,RR_ii,rr_ii;
  const unsigned char* pim;
  int *pnProjRN,*pnProjLN,*numpoR_b,*numpoL_b;

    // check arguments
    if ((buflen==NULL)||(pR==NULL))
      return ERR_GEN_NULLPOINTER;
    // verify radii and calculate max outer radius
    R = *pR;
    if (R<10)
      return ERR_GEN_SIZE_NOT_MATCH;
    if ((r+5>R)||(r<5))
      return ERR_GEN_SIZE_NOT_MATCH;
    if ((xs<2*r)||(ys<2*r))
      return ERR_GEN_NO_DATA;
    if ((xc-r/2<0)||(xc+r/2>=xs)||(yc-r/2<0)||(yc+r/2>=ys))
      return ERR_GEN_INVALID_PARAMETER;
    // verify outer radius - bound by maximum span from center to border
    i = xs-1-xc;
    if (i<xc)
      i = xc;
    if (R>i)
      R = i;
    // verify outer radius - bound by image size
    if (2*R>xs)
      R = xs/2;
    // return corrected outer radius
    *pR = R;
    // calc size
    i = (2+4+2)*(R+1)*sizeof(numpoR[0]);
    if (buf==NULL)
    {
      *buflen = i;
      return ERR_OK;
    }
    if (*buflen<i)
    {
      *buflen = i;
      return ERR_GEN_NOMEMORY;
    }
    // more check arguments
    if ((pnProjRN_b==NULL)||(pnProjLN_b==NULL)||(im==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((((ptr_t)im)&3)||(((ptr_t)buf)&3)||(xs&3)||
        (((ptr_t)pnProjRN_b)&3)||(((ptr_t)pnProjLN_b)&3))
      return ERR_GEN_BAD_ALIGNMENT;
    // allocate
    numpoR   = (int*)buf;
    numpoL   = numpoR+(R+1);
    pnProjRN = numpoL+(R+1);
    pnProjLN = pnProjRN+(R+1);
    numpoR_b = pnProjLN+(R+1);
    numpoL_b = numpoR_b+(R+1);
    // clear junk
    MIA_memset(pnProjRN,0,sizeof(pnProjRN[0])*(1+R));
    MIA_memset(pnProjLN,0,sizeof(pnProjLN[0])*(1+R));
    MIA_memset(numpoR,0,sizeof(numpoR[0])*(1+R));
    MIA_memset(numpoL,0,sizeof(numpoL[0])*(1+R));
    // along y
    if (!nMaxRadOfDence)
      D = 0x10000;
    else
    {
      D = (r*0x10000)/nMaxRadOfDence;
      if (D<0x10000)
        D = 0x10000;
    }
    for (Y=-R*0x10000;Y<=R*0x10000;Y+=D)
    {
      i = (Y+0x8000)/0x10000;
      y = yc+i;
      if ((y<1)||(y>=ys-1))
        continue;
      // left side
      j_out = RR_ii = MIA_INT_sqrt(R*R-i*i);
      if (xc-j_out<1)
        j_out = xc-1;
      j_inn = (i>0)?i:(-i);
      if (j_inn<r)
      {
        rr_ii = MIA_INT_sqrt(r*r-i*i);
        if (j_inn<rr_ii)
          j_inn = rr_ii;
      }
      L = i*i+j_out*j_out;
      currad = MIA_INT_sqrt(L);
      pim = im+y*xs+xc-j_out;
      for (j=-j_out;j<=-j_inn;j++)
      {
        // count number of points in circle
        numpoL[currad]++;
        // straight gradient
        gx =+pim[xs+1]+2*pim[ 1]+pim[-xs+1]
            -pim[xs-1]-2*pim[-1]-pim[-xs-1];
        if (gx<-4)
        {
          gy =+pim[ xs+1]+2*pim[ xs]+pim[ xs-1]
              -pim[-xs+1]-2*pim[-xs]-pim[-xs-1];
          G = gx*gx+gy*gy;
          if (G>2*36)
//          if (G>36)
          {
            // inner product
            S = gx*j+gy*i;
            //
            if ((S>0)&&(((int40)(S*2))*S>((int40)L)*G))
              pnProjLN[currad]++;
          }
        }
        pim++;
        L -= (-j*2-1);
        if (L<currad*currad)
          currad--;
      }
      // right side
      j_out = RR_ii;
      if (xc+j_out>xs-2)
        j_out = xs-2-xc;
      L = i*i+j_out*j_out;
      currad = MIA_INT_sqrt(L);
      pim = im+y*xs+xc+j_out;
      for (j=j_out;j>=j_inn;j--)
      {
        // count number of points in circle
        numpoR[currad]++;
        // straight gradient
        gx =+pim[xs+1]+2*pim[ 1]+pim[-xs+1]
            -pim[xs-1]-2*pim[-1]-pim[-xs-1];
        if (gx>4)
        {
          gy =+pim[ xs+1]+2*pim[ xs]+pim[ xs-1]
              -pim[-xs+1]-2*pim[-xs]-pim[-xs-1];
          G = gx*gx+gy*gy;
          if (G>2*36)
//          if (G>36)
          {
            // inner product
            S = gx*j+gy*i;
            //
            if ((S>0)&&(((int40)(S*2))*S>((int40)L)*G))
              pnProjRN[currad]++;
          }
        }
        pim--;
        L -= (j*2-1);
        if (L<currad*currad)
          currad--;
      }
    }
    // normalize
    for (i=0;i<r;i++)
    {
      numpoR[i] = numpoR[r];
      numpoL[i] = numpoL[r];
      pnProjRN[i] = pnProjRN[r];
      pnProjLN[i] = pnProjLN[r];
    }
    for (i=0;i<=R;i++)
    {
      pnProjRN[i] = (numpoR[i])?((pnProjRN[i]*1024)/numpoR[i]):0;
      pnProjLN[i] = (numpoL[i])?((pnProjLN[i]*1024)/numpoL[i]):0;
    }
    // process histogram: blur
    IPL_HIST_Blur(numpoR_b,numpoR,R+1,BLUR_HIST_WND);
    IPL_HIST_Blur(numpoL_b,numpoL,R+1,BLUR_HIST_WND);
    IPL_HIST_Blur(pnProjRN_b,pnProjRN,R+1,BLUR_HIST_WND);
    IPL_HIST_Blur(pnProjLN_b,pnProjLN,R+1,BLUR_HIST_WND);
    // correct for small presence
    for (i=1;i<=R;i++)
    {
      L = 10000*numpoR_b[i]/(157*(2*BLUR_HIST_WND+1)*i);
      if (L<20)
        pnProjRN_b[i] = 0;
      L = 10000*numpoL_b[i]/(157*(2*BLUR_HIST_WND+1)*i);
      if (L<20)
        pnProjLN_b[i] = 0;
    }
    pnProjRN_b[0] = pnProjRN_b[1];
    pnProjLN_b[0] = pnProjLN_b[1];
    return ERR_OK;
}

// calculate left and right side projection histograms in a concentric ring
MIA_RESULT_CODE IPL_PROJ_FindHoughDonatorProjection5(
  int* pnProjRN_b,          // OUT: circular projection - right side nums
  int* pnProjLN_b,          // OUT: circular projection - left side nums
  const unsigned char* im,  // IN:  image
  int xs,                   // IN:  image size
  int ys,                   // 
  int xc,                   // IN:  ring center
  int yc,                   //
  int r,                    // IN:  inner radius
  int* pR,                  // IN/OUT: outer radius proposed/actually used
  int nMaxRadOfDence,       // IN:  maximum radius for dense processing. 0 - dense always
  void* buf,                // IN:  external buffer
  int* buflen)              // IN/OUT: allocated/used bytes
{
  int i,j,gx,gy,y,G,L,S,*numpoR,*numpoL,currad,R;
  int j_inn,j_out,Y,D,RR_ii,rr_ii;
  const unsigned char* pim;
  int *pnProjRN,*pnProjLN,*numpoR_b,*numpoL_b;

    // check arguments
    if ((buflen==NULL)||(pR==NULL))
      return ERR_GEN_NULLPOINTER;
    // verify radii and calculate max outer radius
    R = *pR;
    if (R<10)
      return ERR_GEN_SIZE_NOT_MATCH;
    if ((r+5>R)||(r<5))
      return ERR_GEN_SIZE_NOT_MATCH;
    if ((xs<2*r)||(ys<2*r))
      return ERR_GEN_NO_DATA;
    if ((xc-r/2<0)||(xc+r/2>=xs)||(yc-r/2<0)||(yc+r/2>=ys))
      return ERR_GEN_INVALID_PARAMETER;
    // verify outer radius - bound by maximum span from center to border
    i = xs-1-xc;
    if (i<xc)
      i = xc;
    if (R>i)
      R = i;
    // verify outer radius - bound by image size
    if (2*R>xs)
      R = xs/2;
    // return corrected outer radius
    *pR = R;
    // calc size
    i = (2+4+2)*(R+1)*sizeof(numpoR[0]);
    if (buf==NULL)
    {
      *buflen = i;
      return ERR_OK;
    }
    if (*buflen<i)
    {
      *buflen = i;
      return ERR_GEN_NOMEMORY;
    }
    // more check arguments
    if ((pnProjRN_b==NULL)||(pnProjLN_b==NULL)||(im==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((((ptr_t)im)&3)||(((ptr_t)buf)&3)||(xs&3)||
        (((ptr_t)pnProjRN_b)&3)||(((ptr_t)pnProjLN_b)&3))
      return ERR_GEN_BAD_ALIGNMENT;
    // allocate
    numpoR   = (int*)buf;
    numpoL   = numpoR+(R+1);
    pnProjRN = numpoL+(R+1);
    pnProjLN = pnProjRN+(R+1);
    numpoR_b = pnProjLN+(R+1);
    numpoL_b = numpoR_b+(R+1);
    // clear junk
    MIA_memset(pnProjRN,0,sizeof(pnProjRN[0])*(1+R));
    MIA_memset(pnProjLN,0,sizeof(pnProjLN[0])*(1+R));
    MIA_memset(numpoR,0,sizeof(numpoR[0])*(1+R));
    MIA_memset(numpoL,0,sizeof(numpoL[0])*(1+R));
    // along y
    if (!nMaxRadOfDence)
      D = 0x10000;
    else
    {
      D = (r*0x10000)/nMaxRadOfDence;
      if (D<0x10000)
        D = 0x10000;
    }
    for (Y=-R*0x10000;Y<=R*0x10000;Y+=D)
    {
      i = (Y+0x8000)/0x10000;
      y = yc+i;
      if ((y<1)||(y>=ys-1))
        continue;
      // left side
      j_out = RR_ii = MIA_INT_sqrt(R*R-i*i);
      if (xc-j_out<1)
        j_out = xc-1;
      j_inn = (i>0)?i:(-i);
      if (j_inn<r)
      {
        rr_ii = MIA_INT_sqrt(r*r-i*i);
        if (j_inn<rr_ii)
          j_inn = rr_ii;
      }
      L = i*i+j_out*j_out;
      currad = MIA_INT_sqrt(L);
      pim = im+y*xs+xc-j_out;
      for (j=-j_out;j<=-j_inn;j++)
      {
        // count number of points in circle
        numpoL[currad]++;
        // straight gradient
        gx =+pim[xs+1]+2*pim[ 1]+pim[-xs+1]
            -pim[xs-1]-2*pim[-1]-pim[-xs-1];
        if (gx<-6)
        {
          gy =+pim[ xs+1]+2*pim[ xs]+pim[ xs-1]
              -pim[-xs+1]-2*pim[-xs]-pim[-xs-1];
          G = gx*gx+gy*gy;
          if (G>2*64)
          {
            // inner product
            S = gx*j+gy*i;
            //
            if ((S>0)&&(((int40)(S*2))*S>((int40)L)*G))
              pnProjLN[currad]++;
          }
        }
        pim++;
        L -= (-j*2-1);
        if (L<currad*currad)
          currad--;
      }
      // right side
      j_out = RR_ii;
      if (xc+j_out>xs-2)
        j_out = xs-2-xc;
      L = i*i+j_out*j_out;
      currad = MIA_INT_sqrt(L);
      pim = im+y*xs+xc+j_out;
      for (j=j_out;j>=j_inn;j--)
      {
        // count number of points in circle
        numpoR[currad]++;
        // straight gradient
        gx =+pim[xs+1]+2*pim[ 1]+pim[-xs+1]
            -pim[xs-1]-2*pim[-1]-pim[-xs-1];
        if (gx>6)
        {
          gy =+pim[ xs+1]+2*pim[ xs]+pim[ xs-1]
              -pim[-xs+1]-2*pim[-xs]-pim[-xs-1];
          G = gx*gx+gy*gy;
          if (G>2*64)
          {
            // inner product
            S = gx*j+gy*i;
            //
            if ((S>0)&&(((int40)(S*2))*S>((int40)L)*G))
              pnProjRN[currad]++;
          }
        }
        pim--;
        L -= (j*2-1);
        if (L<currad*currad)
          currad--;
      }
    }
    // normalize
    for (i=0;i<r;i++)
    {
      numpoR[i] = numpoR[r];
      numpoL[i] = numpoL[r];
      pnProjRN[i] = pnProjRN[r];
      pnProjLN[i] = pnProjLN[r];
    }
    for (i=0;i<=R;i++)
    {
      pnProjRN[i] = (numpoR[i])?((pnProjRN[i]*1024)/numpoR[i]):0;
      pnProjLN[i] = (numpoL[i])?((pnProjLN[i]*1024)/numpoL[i]):0;
    }
    // process histogram: blur
    IPL_HIST_Blur(numpoR_b,numpoR,R+1,BLUR_HIST_WND);
    IPL_HIST_Blur(numpoL_b,numpoL,R+1,BLUR_HIST_WND);
    IPL_HIST_Blur(pnProjRN_b,pnProjRN,R+1,BLUR_HIST_WND);
    IPL_HIST_Blur(pnProjLN_b,pnProjLN,R+1,BLUR_HIST_WND);
    // correct for small presence
    for (i=1;i<=R;i++)
    {
      L = 10000*numpoR_b[i]/(157*(2*BLUR_HIST_WND+1)*i);
      if (L<20)
        pnProjRN_b[i] = 0;
      L = 10000*numpoL_b[i]/(157*(2*BLUR_HIST_WND+1)*i);
      if (L<20)
        pnProjLN_b[i] = 0;
    }
    pnProjRN_b[0] = pnProjRN_b[1];
    pnProjLN_b[0] = pnProjLN_b[1];
    return ERR_OK;
}

// calculate left and right side projection histograms in a concentric ring
MIA_RESULT_CODE IPL_PROJ_FindHoughDonatorProjection5_Mask(
  int* pnProjRN_b,          // OUT: circular projection - right side nums
  int* pnProjLN_b,          // OUT: circular projection - left side nums
  const unsigned char* im,  // IN:  image
  const unsigned char* mask,// IN:  image
  int xs,                   // IN:  image size
  int ys,                   // 
  int xc,                   // IN:  ring center
  int yc,                   //
  int r,                    // IN:  inner radius
  int* pR,                  // IN/OUT: outer radius proposed/actually used
  int nMaxRadOfDence,       // IN:  maximum radius for dense processing. 0 - dense always
  void* buf,                // IN:  external buffer
  int* buflen)              // IN/OUT: allocated/used bytes
{
  int i,j,gx,gy,y,G,L,S,*numpoR,*numpoL,currad,R;
  int j_inn,j_out,Y,D,RR_ii,rr_ii;
  const unsigned char* pim;
  int *pnProjRN,*pnProjLN,*numpoR_b,*numpoL_b;

    // check arguments
    if ((buflen==NULL)||(pR==NULL))
      return ERR_GEN_NULLPOINTER;
    // verify radii and calculate max outer radius
    R = *pR;
    if (R<10)
      return ERR_GEN_SIZE_NOT_MATCH;
    if ((r+5>R)||(r<5))
      return ERR_GEN_SIZE_NOT_MATCH;
    if ((xs<2*r)||(ys<2*r))
      return ERR_GEN_NO_DATA;
    if ((xc-r/2<0)||(xc+r/2>=xs)||(yc-r/2<0)||(yc+r/2>=ys))
      return ERR_GEN_INVALID_PARAMETER;
    // verify outer radius - bound by maximum span from center to border
    i = xs-1-xc;
    if (i<xc)
      i = xc;
    if (R>i)
      R = i;
    // verify outer radius - bound by image size
    if (2*R>xs)
      R = xs/2;
    // return corrected outer radius
    *pR = R;
    // calc size
    i = (2+4+2)*(R+1)*sizeof(numpoR[0]);
    if (buf==NULL)
    {
      *buflen = i;
      return ERR_OK;
    }
    if (*buflen<i)
    {
      *buflen = i;
      return ERR_GEN_NOMEMORY;
    }
    // more check arguments
    if ((pnProjRN_b==NULL)||(pnProjLN_b==NULL)||(im==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((((ptr_t)im)&3)||(((ptr_t)buf)&3)||(xs&3)||
        (((ptr_t)pnProjRN_b)&3)||(((ptr_t)pnProjLN_b)&3))
      return ERR_GEN_BAD_ALIGNMENT;
    // allocate
    numpoR   = (int*)buf;
    numpoL   = numpoR+(R+1);
    pnProjRN = numpoL+(R+1);
    pnProjLN = pnProjRN+(R+1);
    numpoR_b = pnProjLN+(R+1);
    numpoL_b = numpoR_b+(R+1);
    // clear junk
    MIA_memset(pnProjRN,0,sizeof(pnProjRN[0])*(1+R));
    MIA_memset(pnProjLN,0,sizeof(pnProjLN[0])*(1+R));
    MIA_memset(numpoR,0,sizeof(numpoR[0])*(1+R));
    MIA_memset(numpoL,0,sizeof(numpoL[0])*(1+R));
    // along y
    if (!nMaxRadOfDence)
      D = 0x10000;
    else
    {
      D = (r*0x10000)/nMaxRadOfDence;
      if (D<0x10000)
        D = 0x10000;
    }
    for (Y=-R*0x10000;Y<=R*0x10000;Y+=D)
    {
      i = (Y+0x8000)/0x10000;
      y = yc+i;
      if ((y<1)||(y>=ys-1))
        continue;
      // left side
      j_out = RR_ii = MIA_INT_sqrt(R*R-i*i);
      if (xc-j_out<1)
        j_out = xc-1;
      j_inn = (i>0)?i:(-i);
      if (j_inn<r)
      {
        rr_ii = MIA_INT_sqrt(r*r-i*i);
        if (j_inn<rr_ii)
          j_inn = rr_ii;
      }
      L = i*i+j_out*j_out;
      currad = MIA_INT_sqrt(L);
      pim = im+y*xs+xc-j_out;
      for (j=-j_out;j<=-j_inn;j++)
      {
        if (mask[pim-im]==0)
        {
          // count number of points in circle
          numpoL[currad]++;
          // straight gradient
          gx =+pim[xs+1]+2*pim[ 1]+pim[-xs+1]
              -pim[xs-1]-2*pim[-1]-pim[-xs-1];
          if (gx<-6)
          {
            gy =+pim[ xs+1]+2*pim[ xs]+pim[ xs-1]
                -pim[-xs+1]-2*pim[-xs]-pim[-xs-1];
            G = gx*gx+gy*gy;
            if (G>2*64)
            {
              // inner product
              S = gx*j+gy*i;
              //
              if ((S>0)&&(((int40)(S*2))*S>((int40)L)*G))
                pnProjLN[currad]++;
            }
          }
        }
        pim++;
        L -= (-j*2-1);
        if (L<currad*currad)
          currad--;
      }
      // right side
      j_out = RR_ii;
      if (xc+j_out>xs-2)
        j_out = xs-2-xc;
      L = i*i+j_out*j_out;
      currad = MIA_INT_sqrt(L);
      pim = im+y*xs+xc+j_out;
      for (j=j_out;j>=j_inn;j--)
      {
        if (mask[pim-im]==0)
        {
          // count number of points in circle
          numpoR[currad]++;
          // straight gradient
          gx =+pim[xs+1]+2*pim[ 1]+pim[-xs+1]
              -pim[xs-1]-2*pim[-1]-pim[-xs-1];
          if (gx>6)
          {
            gy =+pim[ xs+1]+2*pim[ xs]+pim[ xs-1]
                -pim[-xs+1]-2*pim[-xs]-pim[-xs-1];
            G = gx*gx+gy*gy;
            if (G>2*64)
            {
              // inner product
              S = gx*j+gy*i;
              //
              if ((S>0)&&(((int40)(S*2))*S>((int40)L)*G))
                pnProjRN[currad]++;
            }
          }
        }
        pim--;
        L -= (j*2-1);
        if (L<currad*currad)
          currad--;
      }
    }
    // normalize
    for (i=0;i<r;i++)
    {
      numpoR[i] = numpoR[r];
      numpoL[i] = numpoL[r];
      pnProjRN[i] = pnProjRN[r];
      pnProjLN[i] = pnProjLN[r];
    }
    for (i=0;i<=R;i++)
    {
      pnProjRN[i] = (numpoR[i])?((pnProjRN[i]*1024)/numpoR[i]):0;
      pnProjLN[i] = (numpoL[i])?((pnProjLN[i]*1024)/numpoL[i]):0;
    }
    // process histogram: blur
    IPL_HIST_Blur(numpoR_b,numpoR,R+1,BLUR_HIST_WND);
    IPL_HIST_Blur(numpoL_b,numpoL,R+1,BLUR_HIST_WND);
    IPL_HIST_Blur(pnProjRN_b,pnProjRN,R+1,BLUR_HIST_WND);
    IPL_HIST_Blur(pnProjLN_b,pnProjLN,R+1,BLUR_HIST_WND);
    // correct for small presence
    for (i=1;i<=R;i++)
    {
      L = 10000*numpoR_b[i]/(157*(2*BLUR_HIST_WND+1)*i);
      if (L<20)
        pnProjRN_b[i] = 0;
      L = 10000*numpoL_b[i]/(157*(2*BLUR_HIST_WND+1)*i);
      if (L<20)
        pnProjLN_b[i] = 0;
    }
    pnProjRN_b[0] = pnProjRN_b[1];
    pnProjLN_b[0] = pnProjLN_b[1];
    return ERR_OK;
}

// calculate left side, right side and total circular projection 
// in a concentric ring
MIA_RESULT_CODE IPL_PROJ_FindHoughDonatorProjection3(
  int* pnProjT,             // OUT: circular projection - total
  int* pnProjR,             // OUT: circular projection - right side
  int* pnProjL,             // OUT: circular projection - left side
  const unsigned char* im,  // IN:  image
  int xs,                   // IN:  image size
  int ys,                   // 
  int xc,                   // IN:  ring center
  int yc,                   //
  int r,                    // IN:  inner radius
  int* pR,                  // IN/OUT: outer radius proposed/actually used
  int nMaxRadOfDence,       // IN:  maximum radius for dence processing. 0 - dence always
  void* buf,                // IN:  external buffer
  int* buflen)              // IN/OUT: allocated/used bytes
{
  int i,j,gx,gy,y,G,L,S,*numpoR,*numpoL,*numpoT,currad,R;
  int j_inn,j_out,Y,D,RR_ii,rr_ii;
  const unsigned char* pim;

    // check arguments
    if ((buflen==NULL)||(pR==NULL))
      return ERR_GEN_NULLPOINTER;
    R = *pR;
    if (R<10)
      return ERR_GEN_SIZE_NOT_MATCH;
    // calc size
    i = (R+1)*3*sizeof(numpoT[0]);//+sizeof(short)*xs*ys+sizeof(int)*(xs+ys);
    if (buf==NULL)
    {
      *buflen = i;
      return ERR_OK;
    }
    if (*buflen<i)
    {
      *buflen = i;
      return ERR_GEN_NOMEMORY;
    }
    // more check arguments
    if ((pnProjT==NULL)||(pnProjR==NULL)||(pnProjL==NULL)||(im==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((((ptr_t)im)&3)||(((ptr_t)buf)&3)||(xs&3)||
        (((ptr_t)pnProjT)&3)||(((ptr_t)pnProjR)&3)||(((ptr_t)pnProjL)&3))
      return ERR_GEN_BAD_ALIGNMENT;
    if ((r+5>R)||(r<5))
      return ERR_GEN_SIZE_NOT_MATCH;
    if ((xs<2*r)||(ys<2*r))
      return ERR_GEN_NO_DATA;
    if ((xc-r/2<0)||(xc+r/2>=xs)||(yc-r/2<0)||(yc+r/2>=ys))
      return ERR_GEN_INVALID_PARAMETER;
    // verify outer radius - bound by maximum span from center to border
    i = xs-1-xc;
    if (i<xc)
      i = xc;
    if (R>i)
      R = i;
    // verify outer radius - bound by image size
    if (2*R>ys)
      R = ys/2;
    if (2*R>xs)
      R = xs/2;
    // return corrected outer radius
    *pR = R;
    // allocate
    numpoT = (int*)buf;
    numpoR = numpoT+(R+1);
    numpoL = numpoR+(R+1);
    // clear junk
    MIA_memset(pnProjT,0,sizeof(pnProjT[0])*(1+R));
    MIA_memset(pnProjR,0,sizeof(pnProjR[0])*(1+R));
    MIA_memset(pnProjL,0,sizeof(pnProjL[0])*(1+R));
    MIA_memset(numpoT,0,3*sizeof(numpoT[0])*(1+R));
    // along y
    if (!nMaxRadOfDence)
      D = 0x10000;
    else
    {
      D = (r*0x10000)/nMaxRadOfDence;
      if (D<0x10000)
        D = 0x10000;
    }
    for (Y=-R*0x10000;Y<=R*0x10000;Y+=D)
    {
      i = (Y+0x8000)/0x10000;
      y = yc+i;
      if ((y<1)||(y>=ys-1))
        continue;
      // left side
      j_out = RR_ii = MIA_INT_sqrt(R*R-i*i);
      if (xc-j_out<1)
        j_out = xc-1;
      j_inn = (i>0)?i:(-i);
      if (j_inn<r)
      {
        rr_ii = MIA_INT_sqrt(r*r-i*i);
        if (j_inn<rr_ii)
          j_inn = rr_ii;
      }
      L = i*i+j_out*j_out;
      currad = MIA_INT_sqrt(L);
      pim = im+y*xs+xc-j_out;
      for (j=-j_out;j<=-j_inn;j++)
      {
        // count number of points in circle
        numpoT[currad]++;
        numpoL[currad]++;
        // straight gradient
        gx =+pim[xs+1]+pim[ 1]+pim[-xs+1]
            -pim[xs-1]-pim[-1]-pim[-xs-1];
        if (gx<-4)
        {
          gy =+pim[xs+1]+pim[xs]+pim[xs-1]
              -pim[-xs+1]-pim[-xs]-pim[-xs-1];
          G = gx*gx+gy*gy;
          if (G>36*2)
          {
            // inner product
            S = gx*j+gy*i;
            //
            if ((S>0)&&(((int40)(S*2))*S>((int40)L)*G))
            {
              pnProjT[currad]++;
              pnProjL[currad]++;
            }
          }
        }
        pim++;
        L -= (-j*2-1);
        if (L<currad*currad)
          currad--;
      }
      // right side
      j_out = RR_ii;
      if (xc+j_out>xs-2)
        j_out = xs-2-xc;
      L = i*i+j_out*j_out;
      currad = MIA_INT_sqrt(L);
      pim = im+y*xs+xc+j_out;
      for (j=j_out;j>=j_inn;j--)
      {
        // count number of points in circle
        numpoT[currad]++;
        numpoR[currad]++;
        // straight gradient
        gx =+pim[xs+1]+pim[ 1]+pim[-xs+1]
            -pim[xs-1]-pim[-1]-pim[-xs-1];
        if (gx>4)
        {
          gy =+pim[xs+1]+pim[xs]+pim[xs-1]
              -pim[-xs+1]-pim[-xs]-pim[-xs-1];
          G = gx*gx+gy*gy;
          if (G>36*2)
          {
            // inner product
            S = gx*j+gy*i;
            //
            if ((S>0)&&(((int40)(S*2))*S>((int40)L)*G))
            {
              pnProjT[currad]++;
              pnProjR[currad]++;
            }
          }
        }
        pim--;
        L -= (j*2-1);
        if (L<currad*currad)
          currad--;
      }
    }
    // normalize
    for (i=0;i<=R;i++)
    {
      pnProjT[i] = (pnProjT[i]*1024)/(numpoT[i]+1);
      pnProjR[i] = (pnProjR[i]*1024)/(numpoR[i]+1);
      pnProjL[i] = (pnProjL[i]*1024)/(numpoL[i]+1);
    }
    return ERR_OK;
}

#define BLWIN_I2P 5

// calculate four quadrant projections (for pupil-from-iris)
MIA_RESULT_CODE IPL_PROJ_FindHoughDonatorProjection4_v2(
  int* pnProjR,             // OUT: circular projection - right side
  int* pnProjL,             // OUT: circular projection - left side
  int* pnProjT,             // OUT: circular projection - top side
  int* pnProjB,             // OUT: circular projection - bottom side
  const unsigned char* im,  // IN:  image
  int xs,                   // IN:  image size
  int ys,                   // 
  int xc,                   // IN:  center
  int yc,                   //
  int R,                    // IN:  outer radius
  void* buf,                // IN:  temporary buffer
  int* plen,                // IN/OUT: temporary buffer allocated/used
  const char* nambeg)
{
  MIA_RESULT_CODE res=ERR_OK;
  int *pnProjR_n,*pnProjL_n,*pnProjT_n,*pnProjB_n,*pnProjX,*pnProjX_n;
  int x,y,w,h,cx,cy,gx,gy,cIdx,sz;
  int64 G,L,P,ll;
  uint8 *dbgim=NULL;
//  char nama[FILENAME_MAX];
  double S,Sg,Sgg;
  int T;

    // check arguments
    if (plen==NULL)
      return ERR_GEN_NULLPOINTER;
    if (R<=0)
      return ERR_GEN_NO_DATA;
    // calculate size
    sz = R*6*sizeof(pnProjR_n[0]);
    if (buf==NULL)
    {
      *plen = sz;
      return ERR_OK;
    }
    if (*plen<sz)
      return ERR_GEN_INSUFFICIENT_BUFFER;
    // check more arguments
    if ((pnProjR==NULL)||(pnProjL==NULL)||(pnProjT==NULL)||(pnProjB==NULL)||(im==NULL))
      return ERR_GEN_NULLPOINTER;
    if (((ptr_t)buf)&3)
      return ERR_GEN_BAD_ALIGNMENT;
    if ((xs<2*R)||(ys<2*R))
      return ERR_GEN_BAD_DATA;
    if ((xc-R/10<0)||(xc+R/10>=xs)||(yc-R/10<0)||(yc+R/10>=ys))
      return ERR_GEN_INVALID_PARAMETER; // too close/out of border
    // allocate
    *plen = sz;
    pnProjR_n = (int*)buf;
    pnProjL_n = pnProjR_n+R;
    pnProjT_n = pnProjL_n+R;
    pnProjB_n = pnProjT_n+R;
    pnProjX   = pnProjB_n+R;
    pnProjX_n = pnProjX+R;
    // clear junk
    MIA_memset(pnProjR,0,sizeof(pnProjR[0])*R);
    MIA_memset(pnProjL,0,sizeof(pnProjL[0])*R);
    MIA_memset(pnProjT,0,sizeof(pnProjR[0])*R);
    MIA_memset(pnProjB,0,sizeof(pnProjL[0])*R);
    MIA_memset(pnProjR_n,0,sizeof(pnProjR_n[0])*R*4);
    // obtain dimensions, remember to leave 1-pixel gap for Sobel mask
    x = xc-R;
    y = yc-R;
    w = 2*R+1;
    h = 2*R+1;
    if (x<1)
    {
      w += x-1;
      x = 1;
    }
    if (y<1)
    {
      h += y-1;
      y = 1;
    }
    if (x+w>=xs-1)
      w = xs-x-2;
    if (y+h>=ys-1)
      h = ys-y-2;
    // allocate debug image
    if (nambeg)
    {
      dbgim = (uint8*)malloc(w*h);
      memset(dbgim,0,w*h);
    }
    // calculate gradient threshold
    S = Sg = Sgg = 0.;
    for (cy=y;cy<y+h;cy++)
    {
      for (cx=x;cx<x+w;cx++)
      {
        L = (cx-xc)*(cx-xc)+(cy-yc)*(cy-yc);
        ll = (int)(sqrt((double)L)+.5);
        if ((ll<R/10)||(ll>3*R/4))
          continue;
        gx = im[(cy-1)*xs+(cx+1)]+2*im[(cy  )*xs+(cx+1)]+im[(cy+1)*xs+(cx+1)]
            -im[(cy-1)*xs+(cx-1)]-2*im[(cy  )*xs+(cx-1)]-im[(cy+1)*xs+(cx-1)];
        gy = im[(cy+1)*xs+(cx+1)]+2*im[(cy+1)*xs+(cx  )]+im[(cy+1)*xs+(cx-1)]
            -im[(cy-1)*xs+(cx+1)]-2*im[(cy-1)*xs+(cx  )]-im[(cy-1)*xs+(cx-1)];
        G = gx*gx+gy*gy;
        if (G<100)
          continue;
        P = (cx-xc)*gx+(cy-yc)*gy;
        if (P*P*10<L*G*9)
          continue;
        S++;
        Sg += sqrt((double)G);
        Sgg += G;
      }
    }
    if (S<=1)
    {
      res = ERR_GEN_NO_DATA;
      goto quitti;
    }
    Sg  /= S;
    Sgg /= S;
    //T = (int)(Sg+sqrt(Sgg-Sg*Sg)+.5);
    //T = (int)(Sg+.5);
    T = (int)(Sg-sqrt(Sgg-Sg*Sg)+.5);
    // calculate gradient
    for (cy=y;cy<y+h;cy++)
    {
      for (cx=x;cx<x+w;cx++)
      {
        gx = im[(cy-1)*xs+(cx+1)]+2*im[(cy  )*xs+(cx+1)]+im[(cy+1)*xs+(cx+1)]
            -im[(cy-1)*xs+(cx-1)]-2*im[(cy  )*xs+(cx-1)]-im[(cy+1)*xs+(cx-1)];
        gy = im[(cy+1)*xs+(cx+1)]+2*im[(cy+1)*xs+(cx  )]+im[(cy+1)*xs+(cx-1)]
            -im[(cy-1)*xs+(cx+1)]-2*im[(cy-1)*xs+(cx  )]-im[(cy-1)*xs+(cx-1)];
        L = (cx-xc)*(cx-xc)+(cy-yc)*(cy-yc);
        ll = (int)(sqrt((double)L)+.5);
        if ((cx-xc)+(cy-yc)>0)
        { // right or top
          if ((cx-xc)-(cy-yc)>0)
            pnProjR_n[ll]++;
          else
            pnProjT_n[ll]++;
        }
        else
        { // left or bottom
          if ((cx-xc)-(cy-yc)>0)
            pnProjB_n[ll]++;
          else
            pnProjL_n[ll]++;
        }
        G = gx*gx+gy*gy;
        if (G<100)
          continue;
        P = (cx-xc)*gx+(cy-yc)*gy;
        if (P*P*10<L*G*9)
          continue;
        if (sqrt((double)G)<T)
          continue;
        if ((cx-xc)+(cy-yc)>0)
        { // right or top
          if ((cx-xc)-(cy-yc)>0)
            pnProjR[ll]++;
          else
            pnProjT[ll]++;
        }
        else
        { // left or bottom
          if ((cx-xc)-(cy-yc)>0)
            pnProjB[ll]++;
          else
            pnProjL[ll]++;
        }
        if (dbgim)
          dbgim[(cy-y)*w+(cx-x)] = 255;
      }
    }
    /*if (nambeg)
    {
      sprintf(nama,"%s_1pts.bmp",nambeg);
      DBGL_FILE_SaveUint8Image(dbgim,w,w,h,nama);
    }*/
    // blur & normalize histograms
    IPL_HIST_Blur(pnProjX,pnProjR,R,BLWIN_I2P);
    IPL_HIST_Blur(pnProjX_n,pnProjR_n,R,BLWIN_I2P);
    for (cIdx=0;cIdx<R;cIdx++)
      pnProjR[cIdx] = (pnProjX_n[cIdx])?(1000*pnProjX[cIdx]/pnProjX_n[cIdx]):0;
    for (cIdx=R/10;(pnProjR[cIdx+1]<=pnProjR[cIdx])&&(cIdx<R);cIdx++);
    for (;cIdx>=0;pnProjR[cIdx--]=1);
    for (cIdx=R*3/4;(pnProjR[cIdx-1]<=pnProjR[cIdx])&&(cIdx);cIdx--);
    for (;cIdx<R;pnProjR[cIdx++]=1);
    IPL_HIST_Blur(pnProjX,pnProjL,R,BLWIN_I2P);
    IPL_HIST_Blur(pnProjX_n,pnProjL_n,R,BLWIN_I2P);
    for (cIdx=0;cIdx<R;cIdx++)
      pnProjL[cIdx] = (pnProjX_n[cIdx])?(1000*pnProjX[cIdx]/pnProjX_n[cIdx]):0;
    for (cIdx=R/10;(pnProjL[cIdx+1]<=pnProjL[cIdx])&&(cIdx<R);cIdx++);
    for (;cIdx>=0;pnProjL[cIdx--]=1);
    for (cIdx=R*3/4;(pnProjL[cIdx-1]<=pnProjL[cIdx])&&(cIdx);cIdx--);
    for (;cIdx<R;pnProjL[cIdx++]=1);
    IPL_HIST_Blur(pnProjX,pnProjB,R,BLWIN_I2P);
    IPL_HIST_Blur(pnProjX_n,pnProjB_n,R,BLWIN_I2P);
    for (cIdx=0;cIdx<R;cIdx++)
      pnProjB[cIdx] = (pnProjX_n[cIdx])?(1000*pnProjX[cIdx]/pnProjX_n[cIdx]):0;
    for (cIdx=R/10;(pnProjB[cIdx+1]<=pnProjB[cIdx])&&(cIdx<R);cIdx++);
    for (;cIdx>=0;pnProjB[cIdx--]=1);
    for (cIdx=R*3/4;(pnProjB[cIdx-1]<=pnProjB[cIdx])&&(cIdx);cIdx--);
    for (;cIdx<R;pnProjB[cIdx++]=1);
    IPL_HIST_Blur(pnProjX,pnProjT,R,BLWIN_I2P);
    IPL_HIST_Blur(pnProjX_n,pnProjT_n,R,BLWIN_I2P);
    for (cIdx=0;cIdx<R;cIdx++)
      pnProjT[cIdx] = (pnProjX_n[cIdx])?(1000*pnProjX[cIdx]/pnProjX_n[cIdx]):0;
    for (cIdx=R/10;(pnProjT[cIdx+1]<=pnProjT[cIdx])&&(cIdx<R);cIdx++);
    for (;cIdx>=0;pnProjT[cIdx--]=1);
    for (cIdx=R*3/4;(pnProjT[cIdx-1]<=pnProjT[cIdx])&&(cIdx);cIdx--);
    for (;cIdx<R;pnProjT[cIdx++]=1);
    // free & exit
quitti:;
    if (dbgim)
      free(dbgim);
    return res;
}

// calculate four quadrant projections (for pupil-from-iris)
MIA_RESULT_CODE IPL_PROJ_FindHoughDonatorProjection4(
  int* pnProjR,             // OUT: circular projection - right side
  int* pnProjL,             // OUT: circular projection - left side
  int* pnProjT,             // OUT: circular projection - top side
  int* pnProjB,             // OUT: circular projection - bottom side
  const unsigned char* im,  // IN:  image
  int xs,                   // IN:  image size
  int ys,                   // 
  int xc,                   // IN:  center
  int yc,                   //
  int R,                    // IN:  outer radius
  void* buf,                // IN:  temporary buffer
  int* plen,                // IN/OUT: temporary buffer allocated/used
  const char* nambeg)
{
  MIA_RESULT_CODE res=ERR_OK;
  int *pnProjR_n,*pnProjL_n,*pnProjT_n,*pnProjB_n,*pnProjX,*pnProjX_n;
  int x,y,w,h,cx,cy,gx,gy,cIdx,sz;
  int64 G,L,P,ll;
  uint8 *dbgim=NULL;
//  char nama[FILENAME_MAX];

    // check arguments
    if (plen==NULL)
      return ERR_GEN_NULLPOINTER;
    if (R<=0)
      return ERR_GEN_NO_DATA;
    // calculate size
    sz = R*6*sizeof(pnProjR_n[0]);
    if (buf==NULL)
    {
      *plen = sz;
      return ERR_OK;
    }
    if (*plen<sz)
      return ERR_GEN_INSUFFICIENT_BUFFER;
    // check more arguments
    if ((pnProjR==NULL)||(pnProjL==NULL)||(pnProjT==NULL)||(pnProjB==NULL)||(im==NULL))
      return ERR_GEN_NULLPOINTER;
    if (((ptr_t)buf)&3)
      return ERR_GEN_BAD_ALIGNMENT;
    if ((xs<2*R)||(ys<2*R))
      return ERR_GEN_BAD_DATA;
    if ((xc-R/10<0)||(xc+R/10>=xs)||(yc-R/10<0)||(yc+R/10>=ys))
      return ERR_GEN_INVALID_PARAMETER; // too close/out of border
    // allocate
    *plen = sz;
    pnProjR_n = (int*)buf;
    pnProjL_n = pnProjR_n+R;
    pnProjT_n = pnProjL_n+R;
    pnProjB_n = pnProjT_n+R;
    pnProjX   = pnProjB_n+R;
    pnProjX_n = pnProjX+R;
    // clear junk
    MIA_memset(pnProjR,0,sizeof(pnProjR[0])*R);
    MIA_memset(pnProjL,0,sizeof(pnProjL[0])*R);
    MIA_memset(pnProjT,0,sizeof(pnProjR[0])*R);
    MIA_memset(pnProjB,0,sizeof(pnProjL[0])*R);
    MIA_memset(pnProjR_n,0,sizeof(pnProjR_n[0])*R*4);
    // obtain dimensions, remember to leave 1-pixel gap for Sobel mask
    x = xc-R;
    y = yc-R;
    w = 2*R+1;
    h = 2*R+1;
    if (x<1)
    {
      w += x-1;
      x = 1;
    }
    if (y<1)
    {
      h += y-1;
      y = 1;
    }
    if (x+w>=xs-1)
      w = xs-x-2;
    if (y+h>=ys-1)
      h = ys-y-2;
    // allocate debug image
    if (nambeg)
    {
      dbgim = (uint8*)malloc(w*h);
      memset(dbgim,0,w*h);
    }
    // calculate gradient
    for (cy=y;cy<y+h;cy++)
    {
      for (cx=x;cx<x+w;cx++)
      {
        gx = im[(cy-1)*xs+(cx+1)]+2*im[(cy  )*xs+(cx+1)]+im[(cy+1)*xs+(cx+1)]
            -im[(cy-1)*xs+(cx-1)]-2*im[(cy  )*xs+(cx-1)]-im[(cy+1)*xs+(cx-1)];
        gy = im[(cy+1)*xs+(cx+1)]+2*im[(cy+1)*xs+(cx  )]+im[(cy+1)*xs+(cx-1)]
            -im[(cy-1)*xs+(cx+1)]-2*im[(cy-1)*xs+(cx  )]-im[(cy-1)*xs+(cx-1)];
        L = (cx-xc)*(cx-xc)+(cy-yc)*(cy-yc);
        ll = (int)(sqrt((double)L)+.5);
        if ((cx-xc)+(cy-yc)>0)
        { // right or top
          if ((cx-xc)-(cy-yc)>0)
            pnProjR_n[ll]++;
          else
            pnProjT_n[ll]++;
        }
        else
        { // left or bottom
          if ((cx-xc)-(cy-yc)>0)
            pnProjB_n[ll]++;
          else
            pnProjL_n[ll]++;
        }
        G = gx*gx+gy*gy;
        if (G<200)
          continue;
        P = (cx-xc)*gx+(cy-yc)*gy;
        if (P*P*10<L*G*9)        //if (P*P*12<L*G*11)
          continue;
        if ((cx-xc)+(cy-yc)>0)
        { // right or top
          if ((cx-xc)-(cy-yc)>0)
            pnProjR[ll]++;
          else
            pnProjT[ll]++;
        }
        else
        { // left or bottom
          if ((cx-xc)-(cy-yc)>0)
            pnProjB[ll]++;
          else
            pnProjL[ll]++;
        }
        if (dbgim)
          dbgim[(cy-y)*w+(cx-x)] = 255;
      }
    }
    /*if (nambeg)
    {
      sprintf(nama,"%s_1pts.bmp",nambeg);
      DBGL_FILE_SaveUint8Image(dbgim,w,w,h,nama);
    }*/
    // blur & normalize histograms
    IPL_HIST_Blur(pnProjX,pnProjR,R,BLWIN_I2P);
    IPL_HIST_Blur(pnProjX_n,pnProjR_n,R,BLWIN_I2P);
    for (cIdx=0;cIdx<R;cIdx++)
      pnProjR[cIdx] = (pnProjX_n[cIdx])?(1000*pnProjX[cIdx]/pnProjX_n[cIdx]):0;
    for (cIdx=R/10;(pnProjR[cIdx+1]<=pnProjR[cIdx])&&(cIdx<R);cIdx++);
    for (;cIdx>=0;pnProjR[cIdx--]=1);
    for (cIdx=R*3/4;(pnProjR[cIdx-1]<=pnProjR[cIdx])&&(cIdx);cIdx--);
    for (;cIdx<R;pnProjR[cIdx++]=1);
    IPL_HIST_Blur(pnProjX,pnProjL,R,BLWIN_I2P);
    IPL_HIST_Blur(pnProjX_n,pnProjL_n,R,BLWIN_I2P);
    for (cIdx=0;cIdx<R;cIdx++)
      pnProjL[cIdx] = (pnProjX_n[cIdx])?(1000*pnProjX[cIdx]/pnProjX_n[cIdx]):0;
    for (cIdx=R/10;(pnProjL[cIdx+1]<=pnProjL[cIdx])&&(cIdx<R);cIdx++);
    for (;cIdx>=0;pnProjL[cIdx--]=1);
    for (cIdx=R*3/4;(pnProjL[cIdx-1]<=pnProjL[cIdx])&&(cIdx);cIdx--);
    for (;cIdx<R;pnProjL[cIdx++]=1);
    IPL_HIST_Blur(pnProjX,pnProjB,R,BLWIN_I2P);
    IPL_HIST_Blur(pnProjX_n,pnProjB_n,R,BLWIN_I2P);
    for (cIdx=0;cIdx<R;cIdx++)
      pnProjB[cIdx] = (pnProjX_n[cIdx])?(1000*pnProjX[cIdx]/pnProjX_n[cIdx]):0;
    for (cIdx=R/10;(pnProjB[cIdx+1]<=pnProjB[cIdx])&&(cIdx<R);cIdx++);
    for (;cIdx>=0;pnProjB[cIdx--]=1);
    for (cIdx=R*3/4;(pnProjB[cIdx-1]<=pnProjB[cIdx])&&(cIdx);cIdx--);
    for (;cIdx<R;pnProjB[cIdx++]=1);
    IPL_HIST_Blur(pnProjX,pnProjT,R,BLWIN_I2P);
    IPL_HIST_Blur(pnProjX_n,pnProjT_n,R,BLWIN_I2P);
    for (cIdx=0;cIdx<R;cIdx++)
      pnProjT[cIdx] = (pnProjX_n[cIdx])?(1000*pnProjX[cIdx]/pnProjX_n[cIdx]):0;
    for (cIdx=R/10;(pnProjT[cIdx+1]<=pnProjT[cIdx])&&(cIdx<R);cIdx++);
    for (;cIdx>=0;pnProjT[cIdx--]=1);
    for (cIdx=R*3/4;(pnProjT[cIdx-1]<=pnProjT[cIdx])&&(cIdx);cIdx--);
    for (;cIdx<R;pnProjT[cIdx++]=1);
    // free & exit
    if (dbgim)
      free(dbgim);
    return res;
}

// calculate four quadrant projections (for pupil-from-iris)
MIA_RESULT_CODE IPL_PROJ_FindHoughDonatorProjection4_v1(
  int* pnProjR,             // OUT: circular projection - right side
  int* pnProjL,             // OUT: circular projection - left side
  int* pnProjT,             // OUT: circular projection - top side
  int* pnProjB,             // OUT: circular projection - bottom side
  const unsigned char* im,  // IN:  image
  int xs,                   // IN:  image size
  int ys,                   // 
  int xc,                   // IN:  center
  int yc,                   //
  int R,                    // IN:  outer radius
  void* buf,                // IN:  temporary buffer
  int* plen,                // IN/OUT: temporary buffer allocated/used
  const char* nambeg)
{
  MIA_RESULT_CODE res=ERR_OK;
  int *pnProjR_s,*pnProjL_s,*pnProjT_s,*pnProjB_s,*pnProjX_s,*pnProjX;
  int x,y,w,h,cx,cy,gx,gy,cIdx,sz;
  int64 G,L,P,ll;
  uint8 *dbgim=NULL;
//  char nama[FILENAME_MAX];

    // check arguments
    if (plen==NULL)
      return ERR_GEN_NULLPOINTER;
    if (R<=0)
      return ERR_GEN_NO_DATA;
    // calculate size
    sz = R*6*sizeof(pnProjR[0]);
    if (buf==NULL)
    {
      *plen = sz;
      return ERR_OK;
    }
    if (*plen<sz)
      return ERR_GEN_INSUFFICIENT_BUFFER;
    // check more arguments
    if ((pnProjR==NULL)||(pnProjL==NULL)||(pnProjT==NULL)||(pnProjB==NULL)||(im==NULL))
      return ERR_GEN_NULLPOINTER;
    if (((ptr_t)buf)&3)
      return ERR_GEN_BAD_ALIGNMENT;
    if ((xs<2*R)||(ys<2*R))
      return ERR_GEN_BAD_DATA;
    if ((xc-R/10<0)||(xc+R/10>=xs)||(yc-R/10<0)||(yc+R/10>=ys))
      return ERR_GEN_INVALID_PARAMETER; // too close/out of border
    // allocate
    *plen = sz;
    pnProjR_s = (int*)buf;
    pnProjL_s = pnProjR_s+R;
    pnProjT_s = pnProjL_s+R;
    pnProjB_s = pnProjT_s+R;
    pnProjX_s = pnProjB_s+R;
    pnProjX = pnProjX_s+R;
    // clear junk
    MIA_memset(pnProjR,0,sizeof(pnProjR_s[0])*R);
    MIA_memset(pnProjL,0,sizeof(pnProjL_s[0])*R);
    MIA_memset(pnProjT,0,sizeof(pnProjR_s[0])*R);
    MIA_memset(pnProjB,0,sizeof(pnProjL_s[0])*R);
    MIA_memset(pnProjR_s,0,sizeof(pnProjR[0])*R*4);
    // obtain dimensions, remember to leave 1-pixel gap for Sobel mask
    x = xc-R;
    y = yc-R;
    w = 2*R+1;
    h = 2*R+1;
    if (x<1)
    {
      w += x-1;
      x = 1;
    }
    if (y<1)
    {
      h += y-1;
      y = 1;
    }
    if (x+w>=xs-1)
      w = xs-x-2;
    if (y+h>=ys-1)
      h = ys-y-2;
    // allocate debug image
    if (nambeg)
    {
      dbgim = (uint8*)malloc(w*h);
      memset(dbgim,0,w*h);
    }
    // calculate gradient
    for (cy=y;cy<y+h;cy++)
    {
      for (cx=x;cx<x+w;cx++)
      {
        gx = im[(cy-1)*xs+(cx+1)]+2*im[(cy  )*xs+(cx+1)]+im[(cy+1)*xs+(cx+1)]
            -im[(cy-1)*xs+(cx-1)]-2*im[(cy  )*xs+(cx-1)]-im[(cy+1)*xs+(cx-1)];
        gy = im[(cy+1)*xs+(cx+1)]+2*im[(cy+1)*xs+(cx  )]+im[(cy+1)*xs+(cx-1)]
            -im[(cy-1)*xs+(cx+1)]-2*im[(cy-1)*xs+(cx  )]-im[(cy-1)*xs+(cx-1)];
        L = (cx-xc)*(cx-xc)+(cy-yc)*(cy-yc);
        ll = (int)(sqrt((double)L)+.5);
        if ((cx-xc)+(cy-yc)>0)
        { // right or top
          if ((cx-xc)-(cy-yc)>0)
            pnProjR[ll]++;
          else
            pnProjT[ll]++;
        }
        else
        { // left or bottom
          if ((cx-xc)-(cy-yc)>0)
            pnProjB[ll]++;
          else
            pnProjL[ll]++;
        }
        G = gx*gx+gy*gy;
        if (G<200)
          continue;
        P = (cx-xc)*gx+(cy-yc)*gy;
        if (P*P*10<L*G*9)        //if (P*P*12<L*G*11)
          continue;
        if ((cx-xc)+(cy-yc)>0)
        { // right or top
          if ((cx-xc)-(cy-yc)>0)
            pnProjR_s[ll]++;
          else
            pnProjT_s[ll]++;
        }
        else
        { // left or bottom
          if ((cx-xc)-(cy-yc)>0)
            pnProjB_s[ll]++;
          else
            pnProjL_s[ll]++;
        }
        if (dbgim)
          dbgim[(cy-y)*w+(cx-x)] = 255;
      }
    }
    /*if (nambeg)
    {
      sprintf(nama,"%s_1pts.bmp",nambeg);
      DBGL_FILE_SaveUint8Image(dbgim,w,w,h,nama);
    }*/
    // blur & normalize histograms
    IPL_HIST_Blur(pnProjX_s,pnProjR_s,R,BLWIN_I2P);
    IPL_HIST_Blur(pnProjX,pnProjR,R,BLWIN_I2P);
    for (cIdx=0;cIdx<R;cIdx++)
      pnProjR_s[cIdx] = (pnProjX[cIdx])?(1000*pnProjX_s[cIdx]/pnProjX[cIdx]):0;
    for (cIdx=R/10;(pnProjR_s[cIdx+1]<=pnProjR_s[cIdx])&&(cIdx<R);cIdx++);
    for (;cIdx>=0;pnProjR_s[cIdx--]=1);
    for (cIdx=R*3/4;(pnProjR_s[cIdx-1]<=pnProjR_s[cIdx])&&(cIdx);cIdx--);
    for (;cIdx<R;pnProjR_s[cIdx++]=1);
memcpy(pnProjR,pnProjR_s,R*4);
    IPL_HIST_Blur(pnProjX_s,pnProjL_s,R,BLWIN_I2P);
    IPL_HIST_Blur(pnProjX,pnProjL,R,BLWIN_I2P);
    for (cIdx=0;cIdx<R;cIdx++)
      pnProjL_s[cIdx] = (pnProjX[cIdx])?(1000*pnProjX_s[cIdx]/pnProjX[cIdx]):0;
    for (cIdx=R/10;(pnProjL_s[cIdx+1]<=pnProjL_s[cIdx])&&(cIdx<R);cIdx++);
    for (;cIdx>=0;pnProjL_s[cIdx--]=1);
    for (cIdx=R*3/4;(pnProjL_s[cIdx-1]<=pnProjL_s[cIdx])&&(cIdx);cIdx--);
    for (;cIdx<R;pnProjL_s[cIdx++]=1);
memcpy(pnProjL,pnProjL_s,R*4);
    IPL_HIST_Blur(pnProjX_s,pnProjB_s,R,BLWIN_I2P);
    IPL_HIST_Blur(pnProjX,pnProjB,R,BLWIN_I2P);
    for (cIdx=0;cIdx<R;cIdx++)
      pnProjB_s[cIdx] = (pnProjX[cIdx])?(1000*pnProjX_s[cIdx]/pnProjX[cIdx]):0;
    for (cIdx=R/10;(pnProjB_s[cIdx+1]<=pnProjB_s[cIdx])&&(cIdx<R);cIdx++);
    for (;cIdx>=0;pnProjB_s[cIdx--]=1);
    for (cIdx=R*3/4;(pnProjB_s[cIdx-1]<=pnProjB_s[cIdx])&&(cIdx);cIdx--);
    for (;cIdx<R;pnProjB_s[cIdx++]=1);
memcpy(pnProjB,pnProjB_s,R*4);
    IPL_HIST_Blur(pnProjX_s,pnProjT_s,R,BLWIN_I2P);
    IPL_HIST_Blur(pnProjX,pnProjT,R,BLWIN_I2P);
    for (cIdx=0;cIdx<R;cIdx++)
      pnProjT_s[cIdx] = (pnProjX[cIdx])?(1000*pnProjX_s[cIdx]/pnProjX[cIdx]):0;
    for (cIdx=R/10;(pnProjT_s[cIdx+1]<=pnProjT_s[cIdx])&&(cIdx<R);cIdx++);
    for (;cIdx>=0;pnProjT_s[cIdx--]=1);
    for (cIdx=R*3/4;(pnProjT_s[cIdx-1]<=pnProjT_s[cIdx])&&(cIdx);cIdx--);
    for (;cIdx<R;pnProjT_s[cIdx++]=1);
memcpy(pnProjT,pnProjT_s,R*4);
    // free & exit
    if (dbgim)
      free(dbgim);
    return res;
}
