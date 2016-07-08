/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  28 November 2012                                       */
/*      Revision: 1.0.00                                                 */
/*      Purpose:  diffusion pump                                         */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/
#define __FILENUM__ 114 // __FILENUM__TAG114

#include <stdio.h>
#include <malloc.h>
//#include "dbg.h"
#include "../ownipl.h"

MIA_RESULT_CODE IPL_FILT_AnisotropicDiffusion(
  double* dst,
  const double* src,
  int xs,
  int ys,
  double K,
  double lambda,
  int iterations,
  const char* nambeg)
{
  double tmpdbl;
  double K2Inv = 1/(K*K);
  double* gradX=NULL;
  double* gradY=NULL;
  double* fluxX=NULL;
  double* fluxY=NULL;
  int k,i,j;
//  char nama[FILENAME_MAX];

    memcpy(dst,src,xs*ys*sizeof(dst[0]));
    gradX = (double*)malloc(sizeof(gradX[0])*xs*ys*4);
    gradY = gradX+xs*ys;
    fluxX = gradY+xs*ys;
    fluxY = fluxX+xs*ys;
    memset(gradX,0,sizeof(gradX[0])*xs*ys*4);
    for(k=0;k<iterations;k++)
    {
      // Calculate gradient
      for(i=1;i<ys-1;i++)
        for(j=1;j<xs-1;j++)
        {
          gradX[i*xs+j] = dst[i*xs+j]-dst[i*xs+(j-1)];
          gradY[i*xs+j] = dst[i*xs+j]-dst[(i-1)*xs+j];
        }
/*if (nambeg)
{
  sprintf(nama,"%s_%d_gx.bmp",nambeg,k);
  DBGL_FILE_SaveDoubleImage(gradX,xs,xs,ys,nama);
  sprintf(nama,"%s_%d_gy.bmp",nambeg,k);
  DBGL_FILE_SaveDoubleImage(gradY,xs,xs,ys,nama);
}*/
      // Calculate flux
      for(i=0;i<ys;i++)
        for(j=0;j<xs;j++)
        {
          tmpdbl = gradX[i*xs+j];
          if (tmpdbl<0.)
            tmpdbl = -tmpdbl;
          fluxX[i*xs+j] = gradX[i*xs+j]*exp(-K2Inv*tmpdbl);
          tmpdbl = gradY[i*xs+j];
          if (tmpdbl<0.)
            tmpdbl = -tmpdbl;
          fluxY[i*xs+j] = gradY[i*xs+j]*exp(-K2Inv*tmpdbl);
        }
/*if (nambeg)
{
  sprintf(nama,"%s_%d_fx.bmp",nambeg,k);
  DBGL_FILE_SaveDoubleImage(fluxX,xs,xs,ys,nama);
  sprintf(nama,"%s_%d_fy.bmp",nambeg,k);
  DBGL_FILE_SaveDoubleImage(fluxY,xs,xs,ys,nama);
}*/
      // Update image
      for(i=1;i<ys-1;i++)
      {
        for(j=1;j<xs-1;j++)
        {
          dst[i*xs+j] += lambda*
            (fluxX[i*xs+(j+1)]-fluxX[i*xs+j]+
             fluxY[(i+1)*xs+j]-fluxY[i*xs+j]);
        }
      }
    }
    free(gradX);
    return ERR_OK;
}
