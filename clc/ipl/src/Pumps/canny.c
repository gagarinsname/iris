/*------------------------------------------------------------------------*/
/*                                                                        */
/*      Created:  19 November 2010                                        */
/*      Revision: 1.0.00                                                  */
/*      Purpose:  Canny filter (inspired by Mike Heath's code)            */
/*      Authors:                                                          */
/*        Ivan Matveev / after Mike Heath, refer to                       */
/*  Heath, M., Sarkar, S., Sanocki, T., and Bowyer, K. Comparison of edge */
/*    detectors: a methodology and initial study, Computer Vision and     */
/*    Image Understanding 69 (1), 38-54, January 1998.                    */
/*  Heath, M., Sarkar, S., Sanocki, T. and Bowyer, K.W. A Robust Visual   */
/*    Method for Assessing the Relative Performance of Edge Detection     */
/*    Algorithms, IEEE Transactions on Pattern Analysis and Machine       */
/*    Intelligence 19 (12),  1338-1359, December 1997.                    */
/*                                                                        */
/*------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "stddefs.h"
#include "errorcodes.h"
#include "ipl.h"

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

#define NOEDGE 0 //255
#define POSSIBLE_EDGE 128
#define EDGE 255 //0
#define BOOSTBLURFACTOR 90.0

// recursive (for hysterezis)
staticFunc void ownIPL_follow_edges(
  unsigned char *edgemapptr,
  const short *edgemagptr,
  short lowval,
  int cols)
{
  const short *tempmagptr;
  unsigned char *tempmapptr;
  int i;
  int x[8] = {1,1,0,-1,-1,-1,0,1},
      y[8] = {0,1,1,1,0,-1,-1,-1};

    for(i=0;i<8;i++)
    {
      tempmapptr = edgemapptr - y[i]*cols + x[i];
      tempmagptr = edgemagptr - y[i]*cols + x[i];
      if((*tempmapptr == POSSIBLE_EDGE) && (*tempmagptr > lowval))
      {
        *tempmapptr = (unsigned char) EDGE;
        ownIPL_follow_edges(tempmapptr,tempmagptr, lowval, cols);
      }
    }
}

// hysterezis
staticFunc void ownIPL_apply_hysteresis(
  unsigned char *edge,
  const short *mag,
  const unsigned char *nms,
  int xs,
  int ys,
  double tlow,
  double thigh)
{
  int r, c, pos, numedges, highcount, lowthreshold, highthreshold;
  int hist[32768];
  short maximum_mag=0;

    /****************************************************************************
    * Initialize the edge map to possible edges everywhere the non-maximal
    * suppression suggested there could be an edge except for the border. At
    * the border we say there can not be an edge because it makes the
    * follow_edges_ algorithm more efficient to not worry about tracking an
    * edge off the side of the image.
    ****************************************************************************/
    for(r=0,pos=0;r<ys;r++)
      for(c=0;c<xs;c++,pos++)
        if(nms[pos] == POSSIBLE_EDGE)
          edge[pos] = POSSIBLE_EDGE;
        else
          edge[pos] = NOEDGE;
    for(r=0,pos=0;r<ys;r++,pos+=xs)
    {
      edge[pos] = NOEDGE;
      edge[pos+xs-1] = NOEDGE;
    }
    pos = (ys-1) * xs;
    for(c=0;c<xs;c++,pos++)
    {
      edge[c] = NOEDGE;
      edge[pos] = NOEDGE;
    }
    /****************************************************************************
    * Compute the histogram of the magnitude image. Then use the histogram to
    * compute hysteresis thresholds.
    ****************************************************************************/
    for(r=0;r<32768;r++)
      hist[r] = 0;
    for(r=0,pos=0;r<ys;r++)
      for(c=0;c<xs;c++,pos++)
        if(edge[pos] == POSSIBLE_EDGE)
          hist[mag[pos]]++;
    /****************************************************************************
    * Compute the number of pixels that passed the nonmaximal suppression.
    ****************************************************************************/
    for(r=1,numedges=0;r<32768;r++)
    {
      if(hist[r] != 0)
        maximum_mag = (short)r;
      numedges += hist[r];
    }
    highcount = (int)(numedges * thigh + 0.5);
    /****************************************************************************
    * Compute the high threshold value as the (100 * thigh) percentage point
    * in the magnitude of the gradient histogram of all the pixels that passes
    * non-maximal suppression. Then calculate the low threshold as a fraction
    * of the computed high threshold value. John Canny said in his paper
    * "A Computational Approach to Edge Detection" that "The ratio of the
    * high to low threshold in the implementation is in the range two or three
    * to one." That means that in terms of this implementation, we should
    * choose tlow ~= 0.5 or 0.33333.
    ****************************************************************************/
    r = 1;
    numedges = hist[1];
    while((r<(maximum_mag-1)) && (numedges < highcount))
    {
      r++;
      numedges += hist[r];
    }
    highthreshold = r;
    lowthreshold = (int)(highthreshold * tlow + 0.5);
    /****************************************************************************
    * This loop looks for pixels above the highthreshold to locate edges and
    * then calls follow_edges_ to continue the edge.
    ****************************************************************************/
    for(r=0,pos=0;r<ys;r++)
      for(c=0;c<xs;c++,pos++)
        if((edge[pos] == POSSIBLE_EDGE) && (mag[pos] >= highthreshold))
        {
          edge[pos] = EDGE;
          ownIPL_follow_edges(edge+pos,mag+pos,(short)lowthreshold,xs);
        }
    /****************************************************************************
    * Set all the remaining possible edges to non-edges.
    ****************************************************************************/
    for(r=0,pos=0;r<ys;r++)
      for(c=0;c<xs;c++,pos++)
        if(edge[pos] != EDGE)
          edge[pos] = NOEDGE;
}

// non-maximal suppression
staticFunc void ownIPL_non_max_supp(
  unsigned char *result,
  const short *mag,
  const short *gradx,
  const short *grady,
  int xs,
  int ys)
{
  int rowcount, colcount;
  const short *magrowptr,*magptr;
  const short *gxrowptr,*gxptr;
  const short *gyrowptr,*gyptr;
  int G0,G1,G2,gx,gy,M1,M2;
  unsigned char *resultrowptr, *resultptr;

    memset(result,NOEDGE,xs*ys);
    /****************************************************************************
    * Suppress non-maximum points.
    ****************************************************************************/
    for(rowcount=1,magrowptr=mag+xs+1,gxrowptr=gradx+xs+1,
        gyrowptr=grady+xs+1,resultrowptr=result+xs+1;
        rowcount<ys-2; 
        rowcount++,magrowptr+=xs,gyrowptr+=xs,gxrowptr+=xs,
        resultrowptr+=xs)
    {
      for(colcount=1,magptr=magrowptr,gxptr=gxrowptr,gyptr=gyrowptr,
          resultptr=resultrowptr;
          colcount<xs-2; 
          colcount++,magptr++,gxptr++,gyptr++,resultptr++)
      {
        if ((G0 = *magptr)==0)
          continue;
        gx = *gxptr;
        gy = *gyptr;
        if(gx >= 0)
        {
          if(gy >= 0)
          {
            if (gx >= gy)
            { // 0th octant
              /* 111 */
              /* Left point */
              G1 = *(magptr - 1);
              G2 = *(magptr - xs - 1);
              M1 = (G1-G0)*gx + (G2-G1)*gy;
//M1 = (G1-G0)*gxptr[-1] + (G2-G1)*gyptr[-1];
              /* Right point */
              G1 = *(magptr + 1);
              G2 = *(magptr + xs + 1);
              M2 = (G1-G0)*gx + (G2-G1)*gy;
//M2 = (G1-G0)*gxptr[1] + (G2-G1)*gyptr[1];
            }
            else
            { // 1st octant
              /* 110 */
              /* Left point */
              G1 = *(magptr - xs);
              G2 = *(magptr - xs - 1);
              M1 = (G2-G1)*gx + (G1-G0)*gy;
//M1 = (G2-G1)*gxptr[-xs] + (G1-G0)*gyptr[-xs];
              /* Right point */
              G1 = *(magptr + xs);
              G2 = *(magptr + xs + 1);
              M2 = (G2-G1)*gx + (G1-G0)*gy; 
//M2 = (G2-G1)*gxptr[xs] + (G1-G0)*gyptr[xs];
            }
          }
          else
          { // gx>=0, gy<0
            if (gx >= -gy)
            { // 7th octant
              /* 101 */
              /* Left point */
              G1 = *(magptr - 1);
              G2 = *(magptr + xs - 1);
              M1 = (G1-G0)*gx + (G1-G2)*gy;
//M1 = (G1-G0)*gxptr[-1] + (G1-G2)*gyptr[-1];
              /* Right point */
              G1 = *(magptr + 1);
              G2 = *(magptr - xs + 1);
              M2 = (G1-G0)*gx + (G1-G2)*gy;
//M2 = (G1-G0)*gxptr[1] + (G1-G2)*gyptr[1];
            }
            else
            { // 6th octant
              /* 100 */
              /* Left point */
              G1 = *(magptr + xs);
              G2 = *(magptr + xs - 1);
              M1 = (G2-G1)*gx + (G0-G1)*gy;
//M1 = (G2-G1)*gxptr[xs] + (G0-G1)*gyptr[xs];
              /* Right point */
              G1 = *(magptr - xs);
              G2 = *(magptr - xs + 1);
              M2 = (G2-G1)*gx + (G0-G1)*gy; 
//M2 = (G2-G1)*gxptr[-xs] + (G0-G1)*gyptr[-xs];
            }
          }
        }
        else
        { // gx<0
          if (gy >= 0)
          {
            if (-gx >= gy)
            { // 3rd octant
              /* 011 */
              /* Left point */
              G1 = *(magptr + 1);
              G2 = *(magptr - xs + 1);
              M1 = (G0-G1)*gx + (G2-G1)*gy;
//M1 = (G0-G1)*gxptr[1] + (G2-G1)*gyptr[1];
              /* Right point */
              G1 = *(magptr - 1);
              G2 = *(magptr + xs - 1);
              M2 = (G0-G1)*gx + (G2-G1)*gy;
//M2 = (G0-G1)*gxptr[-1] + (G2-G1)*gyptr[-1];
            }
            else
            { // 2nd octant
              /* 010 */
              /* Left point */
              G1 = *(magptr - xs);
              G2 = *(magptr - xs + 1);
              M1 = (G1-G2)*gx + (G1-G0)*gy;
//M1 = (G1-G2)*gxptr[-xs] + (G1-G0)*gyptr[-xs];
              /* Right point */
              G1 = *(magptr + xs);
              G2 = *(magptr + xs - 1);
              M2 = (G1-G2)*gx + (G1-G0)*gy;
//M2 = (G1-G2)*gxptr[xs] + (G1-G0)*gyptr[xs];
            }
          }
          else
          { // gx<0, gy<0
            if (-gx > -gy)
            { // 4th octant
              /* 001 */
              /* Left point */
              G1 = *(magptr + 1);
              G2 = *(magptr + xs + 1);
              M1 = (G0-G1)*gx + (G1-G2)*gy;
//M1 = (G0-G1)*gxptr[1] + (G1-G2)*gyptr[1];
              /* Right point */
              G1 = *(magptr - 1);
              G2 = *(magptr - xs - 1);
              M2 = (G0-G1)*gx + (G1-G2)*gy;
//M2 = (G0-G1)*gxptr[-1] + (G1-G2)*gyptr[-1];
            }
            else
            { // 5th octant
              /* 000 */
              /* Left point */
              G1 = *(magptr + xs);
              G2 = *(magptr + xs + 1);
              M1 = (G1-G2)*gx + (G0-G1)*gy;
//M1 = (G1-G2)*gxptr[xs] + (G0-G1)*gyptr[xs];
              /* Right point */
              G1 = *(magptr - xs);
              G2 = *(magptr - xs - 1);
              M2 = (G1-G2)*gx + (G0-G1)*gy;
//M2 = (G1-G2)*gxptr[-xs] + (G0-G1)*gyptr[-xs];
            }
          }
        }
        /* Now determine if the current point is a maximum point */
        if ((M1<0)&&(M2<0))
          *resultptr = (unsigned char)POSSIBLE_EDGE;
      }
    }
}

/*******************************************************************************
* Compute the magnitude of the gradient. This is the square root of
* the sum of the squared derivative values.
*******************************************************************************/
staticFunc void ownIPL_magnitude_x_y(
  short *magnitude,
  const short *delta_x,
  const short *delta_y,
  int xs,
  int ys)
{
  int r, c, pos, sq1, sq2;

    for(r=0,pos=0;r<ys;r++)
      for(c=0;c<xs;c++,pos++)
      {
        sq1 = (int)delta_x[pos] * (int)delta_x[pos];
        sq2 = (int)delta_y[pos] * (int)delta_y[pos];
        magnitude[pos] = (short)(0.5 + sqrt((double)sq1 + (double)sq2));
      }
}

/*******************************************************************************
* Compute the first derivative of the image in both the x any y
* directions. The differential filters that are used are:
*                                          -1
*         dx =  -1 0 +1     and       dy =  0
*                                          +1
*******************************************************************************/
staticFunc void ownIPL_derrivative_x_y(
  short *delta_x,
  short *delta_y,
  const short* smoothedim,
  int xs,
  int ys)
{
  int r, c, pos;

    // Compute the x-derivative. Adjust the derivative at the borders to avoid
    // losing pixels.
    for(r=0;r<ys;r++)
    {
      pos = r * xs;
      delta_x[pos] = smoothedim[pos+1] - smoothedim[pos];
      pos++;
      for(c=1;c<(xs-1);c++,pos++)
        delta_x[pos] = smoothedim[pos+1] - smoothedim[pos-1];
      delta_x[pos] = smoothedim[pos] - smoothedim[pos-1];
    }
    // Compute the y-derivative. Adjust the derivative at the borders to avoid
    // losing pixels.
    for(c=0;c<xs;c++)
    {
      pos = c;
      delta_y[pos] = smoothedim[pos+xs] - smoothedim[pos];
      pos += xs;
      for(r=1;r<(ys-1);r++,pos+=xs)
        delta_y[pos] = smoothedim[pos+xs] - smoothedim[pos-xs];
      delta_y[pos] = smoothedim[pos] - smoothedim[pos-xs];
    }
}

// Blur an image with a gaussian filter.
MIA_RESULT_CODE IPL_FILT_GaussianSmooth_uint8(
  short *smoothedim,
  double *tempim,
  const unsigned char *image,
  int xs,
  int ys,
  double sigma)
{
  int r, c, rr, cc,center,i;
  double *kernel,dot,sum,fx;

    if ((smoothedim==NULL)||(tempim==NULL)||(image==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((xs<=0)||(ys<=0)||(sigma<=0.))
      return ERR_GEN_INVALID_PARAMETER;
    // Create a 1-dimensional Haussian smoothing kernel.
    sum = 0.0;
    center = (int)ceil(2.5 * sigma);
    if ((kernel = (double*)calloc(1 + 2*center,sizeof(double))) == NULL)
      return ERR_GEN_NOMEMORY;
    for(i=-center;i<=center;i++)
    {
      fx = (double)(pow(2.71828,-0.5*i*i/(sigma*sigma))/(sigma*sqrt(6.2831853)));
      kernel[i+center] = fx;
      sum += fx;
    }
    for(i=-center;i<=center;i++)
      kernel[i+center] /= sum;
    // Blur in the x - direction.
    for(r=0;r<ys;r++)
      for(c=0;c<xs;c++)
      {
        dot = 0.0;
        sum = 0.0;
        for(cc=(-center);cc<=center;cc++)
          if(((c+cc) >= 0) && ((c+cc) < xs))
          {
            dot += (double)image[r*xs+(c+cc)] * kernel[center+cc];
            sum += kernel[center+cc];
          }
        tempim[r*xs+c] = dot/sum;
      }
    // Blur in the y - direction.
    for(c=0;c<xs;c++)
      for(r=0;r<ys;r++)
      {
        sum = 0.0;
        dot = 0.0;
        for(rr=(-center);rr<=center;rr++)
          if(((r+rr) >= 0) && ((r+rr) < ys))
          {
            dot += tempim[(r+rr)*xs+c] * kernel[center+rr];
            sum += kernel[center+rr];
          }
        smoothedim[r*xs+c] = (short)(dot*BOOSTBLURFACTOR/sum + 0.5);
      }
    free(kernel);
    return ERR_OK;
}

// Blur an image with a gaussian filter
MIA_RESULT_CODE IPL_FILT_GaussianSmooth_ext_uint8(
  short *smoothedim,
  double *tempim,
  const unsigned char *image,
  int xs,
  int ys,
  double sigma_x,
  double sigma_y,
  double _BOOSTBLURFACTOR)
{
  int r, c, rr, cc,center,i;
  double *kernel,dot,sum,fx;

    if ((smoothedim==NULL)||(tempim==NULL)||(image==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((xs<=0)||(ys<=0)||(sigma_x<=0.)||(sigma_y<=0.))
      return ERR_GEN_INVALID_PARAMETER;
    // Blur in the x - direction.
    // Create a 1-dimensional Haussian smoothing kernel.
    sum = 0.0;
    center = (int)ceil(2.5 * sigma_x);
    if ((kernel = (double*)calloc(1 + 2*center,sizeof(double))) == NULL)
      return ERR_GEN_NOMEMORY;
    for(i=-center;i<=center;i++)
    {
      fx = (double)(pow(2.71828,-0.5*i*i/(sigma_x*sigma_x))/(sigma_x*sqrt(6.2831853)));
      kernel[i+center] = fx;
      sum += fx;
    }
    for(i=-center;i<=center;i++)
      kernel[i+center] /= sum;
    // Blur
    for(r=0;r<ys;r++)
      for(c=0;c<xs;c++)
      {
        dot = 0.0;
        sum = 0.0;
        for(cc=(-center);cc<=center;cc++)
          if(((c+cc) >= 0) && ((c+cc) < xs))
          {
            dot += (double)image[r*xs+(c+cc)] * kernel[center+cc];
            sum += kernel[center+cc];
          }
        tempim[r*xs+c] = dot/sum;
      }
    free(kernel);
    // Blur in the y - direction.
    // Create a 1-dimensional Haussian smoothing kernel.
    sum = 0.0;
    center = (int)ceil(2.5 * sigma_y);
    if ((kernel = (double*)calloc(1 + 2*center,sizeof(double))) == NULL)
      return ERR_GEN_NOMEMORY;
    for(i=-center;i<=center;i++)
    {
      fx = (double)(pow(2.71828,-0.5*i*i/(sigma_y*sigma_y))/(sigma_y*sqrt(6.2831853)));
      kernel[i+center] = fx;
      sum += fx;
    }
    for(i=-center;i<=center;i++)
      kernel[i+center] /= sum;
    // Blur
    for(c=0;c<xs;c++)
      for(r=0;r<ys;r++)
      {
        sum = 0.0;
        dot = 0.0;
        for(rr=(-center);rr<=center;rr++)
          if(((r+rr) >= 0) && ((r+rr) < ys))
          {
            dot += tempim[(r+rr)*xs+c] * kernel[center+rr];
            sum += kernel[center+rr];
          }
        smoothedim[r*xs+c] = (short)(dot*_BOOSTBLURFACTOR/sum + 0.5);
      }
    free(kernel);
    return ERR_OK;
}

// Blur an image with a gaussian filter.
MIA_RESULT_CODE IPL_FILT_GaussianSmooth_double(
  double *smoothedim,
  double *tempim,
  const double *image,
  int xs,
  int ys,
  double sigma)
{
  int r, c, rr, cc,center,i;
  double *kernel,dot,sum,fx;

    if ((smoothedim==NULL)||(tempim==NULL)||(image==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((xs<=0)||(ys<=0)||(sigma<=0.))
      return ERR_GEN_INVALID_PARAMETER;
    // Create a 1-dimensional gaussian smoothing kernel.
    sum = 0.0;
    center = (int)ceil(2.5 * sigma);
    if ((kernel = (double*)calloc(1 + 2*center,sizeof(double))) == NULL)
      return ERR_GEN_NOMEMORY;
    for(i=-center;i<=center;i++)
    {
      fx = (double)(pow(2.71828,-0.5*i*i/(sigma*sigma))/(sigma*sqrt(6.2831853)));
      kernel[i+center] = fx;
      sum += fx;
    }
    for(i=-center;i<=center;i++)
      kernel[i+center] /= sum;
    // Blur in the x - direction.
    for(r=0;r<ys;r++)
      for(c=0;c<xs;c++)
      {
        dot = 0.0;
        sum = 0.0;
        for(cc=(-center);cc<=center;cc++)
          if(((c+cc) >= 0) && ((c+cc) < xs))
          {
            dot += image[r*xs+(c+cc)] * kernel[center+cc];
            sum += kernel[center+cc];
          }
        tempim[r*xs+c] = dot/sum;
      }
    // Blur in the y - direction.
    for(c=0;c<xs;c++)
      for(r=0;r<ys;r++)
      {
        sum = 0.0;
        dot = 0.0;
        for(rr=(-center);rr<=center;rr++)
          if(((r+rr) >= 0) && ((r+rr) < ys))
          {
            dot += tempim[(r+rr)*xs+c] * kernel[center+rr];
            sum += kernel[center+rr];
          }
        smoothedim[r*xs+c] = dot*BOOSTBLURFACTOR/sum;
      }
    free(kernel);
    return ERR_OK;
}

// canny edge detection
//   1) Convolve the image with a separable gaussian filter.
//   2) Take the dx and dy the first derivatives using [-1,0,1] and [1,0,-1]'.
//   3) Compute the magnitude: sqrt(dx*dx+dy*dy).
//   4) Perform non-maximal suppression.
//   5) Perform hysteresis.
//   sigma = The standard deviation of the gaussian smoothing filter.
//   tlow  = Specifies the low value to use in hysteresis. This is a 
//           fraction (0-1) of the computed high threshold edge strength value.
//   thigh = Specifies the high value to use in hysteresis. This fraction (0-1)
//           specifies the percentage point in a histogram of the gradient of
//           the magnitude. Magnitude values of zero are not counted in the
//           histogram.
MIA_RESULT_CODE IPL_FILT_Canny(
  unsigned char *edge,        // OUT: canny mask
  const unsigned char *image, // IN:  source image
  short ** ppgx,              // OUT: map of horizontal gradients
  short ** ppgy,              // OUT: map of vertical gradients
  int xs,                     // IN:  image size
  int ys,                     // 
  double sigma,               // IN:  haussian sigma
  double tlow,                // IN:  low hysteresis threshold (turning off)
  double thigh,               // IN:  high hysteresis threshold (turning on)
  void* pvBuf,                // IN:  temporary buffer
  int* pnLen,                 // IN/OUT: bytes allocated/used
  const char *nambeg)         // IN:  set this to NULL
{
  unsigned char *nms;
  short *smoothedim,*delta_x,*delta_y,*magnitude;
  double *tempim;
//  void* buf=NULL;
  int sz;

    if ((xs<=0)||(ys<=0))
      return ERR_GEN_NO_DATA;
    if (pnLen==NULL)
      return ERR_GEN_NULLPOINTER;
    sz = xs*ys*sizeof(short)+   // smoothed
         xs*ys*sizeof(double)+   // temp double
         xs*ys*sizeof(short)+   // delta_x
         xs*ys*sizeof(short)+   // delta_y
         xs*ys*sizeof(short)+   // magnitude
         xs*ys;                 // nms
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
    smoothedim = (short*)pvBuf;
    delta_x = smoothedim+xs*ys;
    delta_y = delta_x+xs*ys;
    magnitude = delta_y+xs*ys;
    nms = (unsigned char*)(magnitude+xs*ys);
    tempim = (double*)(nms+xs*ys);
    if (ppgx)
      *ppgx = delta_x;
    if (ppgy)
      *ppgy = delta_y;
    // Perform Haussian smoothing
    IPL_FILT_GaussianSmooth_uint8(smoothedim,tempim,image,xs,ys,sigma);
//    gaussian_smooth_int(smoothedim,(short*)tempim,image,xs,ys,sigma);
//IPL_FILT_GaussianSmooth_uint8(delta_x,tempim,image,xs,ys,sigma);
//for (i=xs*ys-1;i>=0;i--)
//  delta_y[i] = smoothedim[i]-delta_x[i];
//sprintf(nam,"%s_1diff.bmp",nambeg);
//DBGL_FILE_SaveInt16Image(delta_y,xs,ys,nam);
//    IPL_FILT_HaussBlur3x3(nms,image,xs,ys);
//    IPL_FILT_CalculateMean(nms,xs,image,xs,xs,ys,0,0,xs,ys,13,13);
//    for (i=xs*ys-1;i>=0;i--)
  //    smoothedim[i] = nms[i];
//sprintf(nam,"%s_1hauss.bmp",nambeg);
//DBGL_FILE_SaveInt16Image(smoothedim,xs,ys,nam);
    // Compute the first derivative in the x and y directions.
    ownIPL_derrivative_x_y(delta_x,delta_y,smoothedim,xs,ys);
//DBGL_FILE_SaveInt16Image(delta_x,xs,ys,"2a-dx.bmp");
//DBGL_FILE_SaveInt16Image(delta_y,xs,ys,"2b-dy.bmp");
    // Compute the magnitude of the gradient.
    ownIPL_magnitude_x_y(magnitude,delta_x,delta_y,xs,ys);
//DBGL_FILE_SaveInt16Image(magnitude,xs,ys,"3-mag.bmp");
    // Perform non-maximal suppression.
    ownIPL_non_max_supp(nms,magnitude,delta_x,delta_y,xs,ys);
//sprintf(nam,"%s_3nms.bmp",nambeg);
//DBGL_FILE_SaveUint8Image(nms,xs,xs,ys,nam);
    // Use hysteresis to mark the edge pixels.
    ownIPL_apply_hysteresis(edge,magnitude,nms,xs,ys,tlow,thigh);
    return ERR_OK;
}
