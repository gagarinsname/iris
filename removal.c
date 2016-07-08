#define __FILENUM__ 999 // __FILENUM__TAG43

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "ipl.h"
#include "bdl.h"
#include "ffl.h"

static int IP_IVIR_RemoveFlashes(
                 unsigned char* im,    // IN:  input image
                 int xs,                     // IN:  image size
                 int ys   );
static int IP_IVIR_RemoveFlashesEx(
                 unsigned char* im,    // IN:  input image
                 int xs,                     // IN:  image size
                 int ys, 
				 unsigned char* pFlashMask   );
static int IP_IVIR_DetectFlashes(
              unsigned char* im,    // IN:  input image
              int xs,                     // IN:  image size
              int ys,
              unsigned char* mask);    // OUT:  mask

#ifndef PI
#define PI         (3.14159265358979323846)
#endif

#ifndef __min
#define __min(a,b) (a>b?b:a)
#endif


void IP_IVIR_PMP_Erode3x3(
              unsigned char* dst,
              const unsigned char* src,
              int xs,
              int ys);
void IP_IVIR_PMP_Erode3x3Cross(
                 unsigned char* dst,
                 const unsigned char* src,
                 int xs,
                 int ys);
void IP_IVIR_PMP_Dilate3x3(
               unsigned char* dst,
               const unsigned char* src,
               int xs,
               int ys);
void IP_IVIR_PMP_Dilate3x3Cross(
                unsigned char* dst,
                const unsigned char* src,
                int xs,
                int ys);

void IP_IVIR_PMP_FillBorders1B(
                 unsigned char* im,
                 int xs,
                 int ys);

/*
#define SAVEMATRIX(matrix, sizeX, sizeY) {\
  int i, j;\
  FILE* out = fopen("c:\\\\" #matrix "Libdll.txt", "wt");\
  fprintf( out, "%d %d", (int)sizeY, (int)sizeX); fprintf( out, "\n");\
  for(j=0; j<sizeY; j++){\
  for(i=0; i<sizeX; i++)\
  fprintf( out, "%9.4f ", (double)((matrix)[i+j*sizeX]) );\
  fprintf( out, "\n");\
  }\
  fclose( out );\
}
*/
#define SAVEMATRIX(matrix, sizeX, sizeY) {\
	FILE* out = fopen("c:\\\\" #matrix "Libdll.txt", "wt");\
	fprintf( out, "%d %d", (int)sizeY, (int)sizeX); fprintf( out, "\n");\
	for(j=0; j<sizeY; j++){\
		for(i=0; i<sizeX; i++)\
			fprintf( out, "%9.4f ", (double)((matrix)[i+j*sizeX]) );\
		fprintf( out, "\n");\
	}\
	fclose( out );\
}

MIA_RESULT_CODE /*FLA_RemoveFlash*/BDL_Flash_RemoveFlash(unsigned char* pImage, int nWidth, int nHeight)
{
  IP_IVIR_RemoveFlashes(pImage, nWidth, nHeight);
  return ERR_OK;
}


MIA_RESULT_CODE /*FLA_RemoveFlash*/BDL_Flash_RemoveFlashEx(unsigned char* pImage, int nWidth, int nHeight, unsigned char* pFlashMask)
{
  IP_IVIR_RemoveFlashesEx(pImage, nWidth, nHeight, pFlashMask);
  return ERR_OK;
}

MIA_RESULT_CODE BDL_Flash_DetectFlash(unsigned char* pImage, int nWidth, int nHeight, unsigned char* pFlashMask)
{
  IP_IVIR_DetectFlashes(pImage, nWidth, nHeight, pFlashMask);
  return ERR_OK;
}

// remove flashes by harmonic filling
int IP_IVIR_RemoveFlashesEx(
              unsigned char* im,    // IN:  input image
              int xs,                     // IN:  image size
              int ys, 
			  unsigned char* pFlashMask   )  // mask of flash to remove
{
  int ir, ic, numIter;
//  int threshDiff2;
  double fi=1.65;		// optimization factor
  int sv_f = 0;      // save or not results
  int i,j;


  //reserve memory for processing of data
//  int histVector[256];                   // histogram array for threshhold calculation
  unsigned char* imPnt = &im[0];

  // optimized reservation
  unsigned char* mask_flashes0;
  unsigned char* mask_flashes;
  unsigned char* Im_eroded_1;
  unsigned char* Im_eroded_2;
  unsigned char* Im_dilated_1;
  unsigned char* Im_dilated_2;
  unsigned char* Im_diff;
  unsigned char* imgRes2;  // for reduced image in 2 times
  double* Img_Iter1; // for iteration  during solution poison equation
  double* Img_Iter2; // for itaration durin solution poison equation
  unsigned char* mask_flashesTmp;

//  unsigned char* mask_flashes_Purified;
  unsigned char* img_Final;
  unsigned char* mask_flashes_Final;

  mask_flashes0 = (unsigned char*) malloc(xs/2*ys/2*8 + xs/2*ys/2*3*8);
  if (!mask_flashes0)
  {
    return -1;
  }
  mask_flashes = (unsigned char*) (mask_flashes0+xs/2*ys/2);
  Im_eroded_1 = (unsigned char*) (mask_flashes+xs/2*ys/2);
  Im_eroded_2 = (unsigned char*) (Im_eroded_1+xs/2*ys/2);
  Im_dilated_1 = (unsigned char*) (Im_eroded_2+xs/2*ys/2);
  Im_dilated_2 = (unsigned char*) (Im_dilated_1+xs/2*ys/2);
  Im_diff = (unsigned char*) (Im_dilated_2+xs/2*ys/2);
  imgRes2 = (unsigned char*) (Im_diff+xs/2*ys/2 );  // for reduced image in 2 times
  Img_Iter1 = (double*) (imgRes2+xs/2*ys/2 ); // for iteration  during solution poison equation
  Img_Iter2 = (double*) (Img_Iter1+xs/2*ys/2 ); // for itaration durin solution poison equation
  mask_flashesTmp = (unsigned char*) (Img_Iter2+xs/2*ys/2);

  // reduce size 2 times
  for( ir=0; ir<ys/2; ir++)
  {
    for( ic=0; ic<xs/2; ic++)
    {
      imgRes2[ir*xs/2+ic] = (unsigned char)( 0.25*( im[ir*2*xs+ic*2] + im[ir*2*xs+ic*2+1] + 
        im[(ir*2+1)*xs+ic*2] + im[(ir*2+1)*xs+ic*2+1] + 2 ) );
      //imgRes2[ir*xs/2+ic] = __min(
      //	__min( im[ir*2*xs+ic*2], im[ir*2*xs+ic*2+1]),
      //	__min(im[(ir*2+1)*xs+ic*2], im[(ir*2+1)*xs+ic*2+1] ) 
      //	) ;
      Img_Iter2[ir*xs/2+ic] = imgRes2[ir*xs/2+ic];

	  mask_flashes[ir*xs/2+ic] = pFlashMask[ir*2*xs+ic*2];
    }
  }
  //if(sv_f) SAVEMATRIX(im,xs,ys);
  if(sv_f) SAVEMATRIX(imgRes2,xs/2,ys/2);
  if(sv_f) SAVEMATRIX(mask_flashes,xs/2,ys/2);



  // solving Poisson's Equation
  // optimization factor for sucsesive relaxation
  //double fi=2/(1+(1-(1-2*sin(3.1415926/2/15).^2).^2).^0.5); 
  for( numIter=0; numIter<10; numIter++)
  { 

    for( ir=1; ir<ys/2-1; ir++)
    {
      for( ic=1; ic<xs/2-1; ic++)
      {

//        if( 1 )// use acceleration    
          if( mask_flashes[ir*xs/2+ic]!=0 )
          {
            Img_Iter1[ir*xs/2+ic] = 0.25*( Img_Iter2[(ir-1)*xs/2+ic] + Img_Iter2[ir*xs/2+ic-1] +
              Img_Iter2[(ir+1)*xs/2+ic] + Img_Iter2[ir*xs/2+ic+1] );
            Img_Iter2[ir*xs/2+ic] = fi* Img_Iter1[ir*xs/2+ic] + (1-fi)*Img_Iter2[ir*xs/2+ic];

          }
//
  //        else
    //        if( mask_flashes[ir*xs/2+ic]!=0 )
      //        Img_Iter2[ir*xs/2+ic] = 0.25*( Img_Iter2[(ir-1)*xs/2+ic] + Img_Iter2[ir*xs/2+ic-1] +
        //      Img_Iter2[(ir+1)*xs/2+ic] + Img_Iter2[ir*xs/2+ic+1] );

      }
    }

  }


      // copy result to image (whith resizing) 
      for( ir=1; ir<ys/2-1; ir++)
      {
        for( ic=1; ic<xs/2-1; ic++)
        {
          if( mask_flashes[ir*xs/2+ic]!=0 )
          {
            double tmpd = Img_Iter2[ir*xs/2+ic];
            if(tmpd<0) tmpd = 0;
            if(tmpd>255) tmpd = 255;
            //if(imPnt[ir*2*xs+ic*2] > tmpd)
            imPnt[ir*2*xs+ic*2] = (unsigned char)tmpd;
            //if(imPnt[ir*2*xs+ic*2+1] > tmpd)
            imPnt[ir*2*xs+ic*2+1] = (unsigned char)tmpd;
            //if(imPnt[(ir*2+1)*xs+ic*2] > tmpd)
            imPnt[(ir*2+1)*xs+ic*2] = (unsigned char)tmpd;
            //if(imPnt[(ir*2+1)*xs+ic*2+1] > tmpd)
            imPnt[(ir*2+1)*xs+ic*2+1] = (unsigned char)tmpd;


          }

//          if(0)
  //        {
    //        imPnt[ir*2*xs+ic*2] = mask_flashes[ir*xs/2+ic];
      //      imPnt[ir*2*xs+ic*2+1] = mask_flashes[ir*xs/2+ic];
        //    imPnt[(ir*2+1)*xs+ic*2] = mask_flashes[ir*xs/2+ic];
          //  imPnt[(ir*2+1)*xs+ic*2+1] = mask_flashes[ir*xs/2+ic];			
          //}
        }
      }

      img_Final = imPnt;
      mask_flashes_Final = mask_flashes;
      if(sv_f) SAVEMATRIX(img_Final,xs,ys);
      if(sv_f) SAVEMATRIX(mask_flashes_Final,xs/2,ys/2);
      // free all memory for processing
      free(mask_flashes0);
      // return DLLT_OK;
  return 0;
}

// remove flashes by harmonic filling
int IP_IVIR_RemoveFlashes(
              unsigned char* im,    // IN:  input image
              int xs,                     // IN:  image size
              int ys   )
{
  int ir, ic, numIter, tmpint;
  int iAccum, threshDiff1, threshDiff2, threshDiff12;
  double fi=1.65;		// optimization factor
  int sv_f = 0;      // save or not results
  int i,j;


  //reserve memory for processing of data
  int histVector[256];                   // histogram array for threshhold calculation
  unsigned char* imPnt = &im[0];

  // optimized reservation
  unsigned char* mask_flashes0;
  unsigned char* mask_flashes;
  unsigned char* Im_eroded_1;
  unsigned char* Im_eroded_2;
  unsigned char* Im_dilated_1;
  unsigned char* Im_dilated_2;
  unsigned char* Im_diff;
  unsigned char* imgRes2;  // for reduced image in 2 times
  double* Img_Iter1; // for iteration  during solution poison equation
  double* Img_Iter2; // for itaration durin solution poison equation
  unsigned char* mask_flashesTmp;

  unsigned char* mask_flashes_Purified;
  unsigned char* img_FirstIter;
  unsigned char* mask_flashes_ExcBr;
  unsigned char* img_Final;
  unsigned char* mask_flashes_Final;

  mask_flashes0 = (unsigned char*) malloc(xs/2*ys/2*8 + xs/2*ys/2*3*8);
  if (!mask_flashes0)
  {
    return -1;
  }
  mask_flashes = (unsigned char*) (mask_flashes0+xs/2*ys/2);
  Im_eroded_1 = (unsigned char*) (mask_flashes+xs/2*ys/2);
  Im_eroded_2 = (unsigned char*) (Im_eroded_1+xs/2*ys/2);
  Im_dilated_1 = (unsigned char*) (Im_eroded_2+xs/2*ys/2);
  Im_dilated_2 = (unsigned char*) (Im_dilated_1+xs/2*ys/2);
  Im_diff = (unsigned char*) (Im_dilated_2+xs/2*ys/2);
  imgRes2 = (unsigned char*) (Im_diff+xs/2*ys/2 );  // for reduced image in 2 times
  Img_Iter1 = (double*) (imgRes2+xs/2*ys/2 ); // for iteration  during solution poison equation
  Img_Iter2 = (double*) (Img_Iter1+xs/2*ys/2 ); // for itaration durin solution poison equation
  mask_flashesTmp = (unsigned char*) (Img_Iter2+xs/2*ys/2);

  // reduce size 2 times
  for( ir=0; ir<ys/2; ir++)
  {
    for( ic=0; ic<xs/2; ic++)
    {
      imgRes2[ir*xs/2+ic] = (unsigned char)( 0.25*( im[ir*2*xs+ic*2] + im[ir*2*xs+ic*2+1] + 
        im[(ir*2+1)*xs+ic*2] + im[(ir*2+1)*xs+ic*2+1] + 2 ) );
      //imgRes2[ir*xs/2+ic] = __min(
      //	__min( im[ir*2*xs+ic*2], im[ir*2*xs+ic*2+1]),
      //	__min(im[(ir*2+1)*xs+ic*2], im[(ir*2+1)*xs+ic*2+1] ) 
      //	) ;
      Img_Iter2[ir*xs/2+ic] = imgRes2[ir*xs/2+ic];
    }
  }
  //if(sv_f) SAVEMATRIX(im,xs,ys);
  if(sv_f) SAVEMATRIX(imgRes2,xs/2,ys/2);


  // erosion-dilation

  // erosion 5 times
  memcpy(Im_eroded_1, imgRes2, xs/2*ys/2);
  for(ic=0; ic<2; ic++)
  {
    IP_IVIR_PMP_Erode3x3(Im_eroded_2,Im_eroded_1,xs/2,ys/2);
    IP_IVIR_PMP_FillBorders1B(  Im_eroded_2,xs/2,ys/2);     // removing borders gap
    memcpy(Im_eroded_1, Im_eroded_2, xs/2*ys/2);
    IP_IVIR_PMP_Erode3x3Cross(Im_eroded_2,Im_eroded_1,xs/2,ys/2);
    IP_IVIR_PMP_FillBorders1B(  Im_eroded_2,xs/2,ys/2);     // removing borders gap
    memcpy(Im_eroded_1, Im_eroded_2, xs/2*ys/2);
  }
  if(sv_f) SAVEMATRIX(Im_eroded_1,xs/2,ys/2);

  // dilation 5 times
  memcpy(Im_dilated_1, Im_eroded_1, xs/2*ys/2);
  for(ic=0; ic<2; ic++)
  {
    IP_IVIR_PMP_Dilate3x3Cross(Im_dilated_2,Im_dilated_1,xs/2,ys/2);
    memcpy(Im_dilated_1, Im_dilated_2, xs/2*ys/2);
    IP_IVIR_PMP_Dilate3x3(Im_dilated_2,Im_dilated_1,xs/2,ys/2);
    memcpy(Im_dilated_1, Im_dilated_2, xs/2*ys/2);
  }
  if(sv_f) SAVEMATRIX(Im_dilated_1,xs/2,ys/2);


  // calculate opening: difference of image and opening
  // and fill histogramm
  memset(histVector, 0, 256*sizeof(int) );
  for( ir=0; ir<ys/2*xs/2; ir++)
  {
    tmpint = imgRes2[ir]-Im_eroded_2[ir]; 
    if(tmpint<0) tmpint = 0; 
    if(tmpint>255) tmpint = 255; 

    Im_diff[ir] = (unsigned char)tmpint;
    //if()
    histVector[Im_diff[ir]]++; 
  }
  if(sv_f) SAVEMATRIX(Im_diff,xs/2,ys/2);

  // calculate threshhold as 
  //threshDiff1 = prctile(imgDiff(:), 100.*(1.-50/prod(size(imgDiff)) )  );
  //threshDiff2 = prctile(imgDiff(:), 99.5  );
  //threshDiff2 = threshDiff1-10;
  //threshDiff12 = ( threshDiff1+2*threshDiff2)/3;

  iAccum = 0;
  for( ir=255; ir>0; ir--)
  {
    iAccum+=histVector[ir];
    if(iAccum>50) // 50 pixels from top
      break;
  }
  threshDiff1 = ir;

  iAccum = 0;
  for( ir=255; ir>0; ir--)
  {
    iAccum+=histVector[ir];
    if(iAccum>ys/2*xs/2/200) // 99.5%
      break;
  }
  threshDiff2 = ir-10;
  threshDiff12 = ( threshDiff1+2*threshDiff2)/3;

  // use threshold for masking
  for( ir=0; ir<ys/2*xs/2; ir++)
  {
    mask_flashes0[ir] = 0;
    if( Im_diff[ir]>threshDiff12 )
      mask_flashes0[ir] = 255;
  }
  //if(sv_f) ; VCL_BMP_Save(mask_flashes0,xs/2,ys/2,xs/2,"c:\\Im_t_masked_0.bmp",USUAL_BMP);

  // dilate mask
  IP_IVIR_PMP_Dilate3x3Cross(mask_flashes, mask_flashes0,xs/2,ys/2);
  memcpy(mask_flashes0, mask_flashes, xs/2*ys/2);
  IP_IVIR_PMP_Dilate3x3(mask_flashes, mask_flashes0,xs/2,ys/2);
  memcpy(mask_flashes0, mask_flashes, xs/2*ys/2);
  IP_IVIR_PMP_Dilate3x3Cross(mask_flashes, mask_flashes0,xs/2,ys/2);

  if(sv_f) SAVEMATRIX(mask_flashes0,xs/2,ys/2);


  // tune for brightness
  for( ir=0; ir<ys/2*xs/2; ir++)
  {
    mask_flashesTmp[ir] = 0;
    if( imgRes2[ir]>200 )
      mask_flashesTmp[ir] = 255;
  }

  // dilate mask
  IP_IVIR_PMP_Dilate3x3(mask_flashes0, mask_flashesTmp,xs/2,ys/2);
  memcpy(mask_flashesTmp, mask_flashes0, xs/2*ys/2);
  IP_IVIR_PMP_Dilate3x3Cross(mask_flashes0, mask_flashesTmp,xs/2,ys/2);
  memcpy(mask_flashesTmp, mask_flashes0, xs/2*ys/2);
  IP_IVIR_PMP_Dilate3x3(mask_flashes0, mask_flashesTmp,xs/2,ys/2);
  memcpy(mask_flashesTmp, mask_flashes0, xs/2*ys/2);
  IP_IVIR_PMP_Dilate3x3Cross(mask_flashes0, mask_flashesTmp,xs/2,ys/2);
  memcpy(mask_flashesTmp, mask_flashes0, xs/2*ys/2);
  IP_IVIR_PMP_Dilate3x3(mask_flashes0, mask_flashesTmp,xs/2,ys/2);
  memcpy(mask_flashesTmp, mask_flashes0, xs/2*ys/2);
  IP_IVIR_PMP_Dilate3x3Cross(mask_flashes0, mask_flashesTmp,xs/2,ys/2);

  for( ir=0; ir<ys/2*xs/2; ir++)
    mask_flashes[ir] = __min( mask_flashes0[ir], mask_flashes[ir] );

  if(sv_f) SAVEMATRIX(mask_flashes,xs/2,ys/2);


  // correct mask to exclude minimal brightness
  for( ir=3; ir<ys/2-3; ir++)
  {
    for( ic=3; ic<xs/2-3; ic++)
    {
      int isVertMaskBr;
      int isHorzMaskBr;
      int isBr;
      double minBr;
      double maxBr;
      int ir1;
      int ic1;
      int indx = ir*xs/2+ic;
      mask_flashesTmp[indx] = 0;

      isVertMaskBr = mask_flashes[indx]+mask_flashes[indx-xs/2]+mask_flashes[indx+xs/2];
      isHorzMaskBr = mask_flashes[indx]+mask_flashes[indx-1]+mask_flashes[indx+1];
      isBr = 0;
      isVertMaskBr /= 255;
      isHorzMaskBr /= 255;
      if(isVertMaskBr==0 || isVertMaskBr==3)
        continue;
      if(isHorzMaskBr==0 || isHorzMaskBr==3)
        continue;

      minBr = 10000;
      maxBr = 0;
      for(ir1=-2;ir1<=2;ir1++)
      {
        for(ic1=-2;ic1<=2;ic1++)
        {
          if(minBr>imgRes2[indx+ic1+ir1*xs/2])
            minBr=imgRes2[indx+ic1+ir1*xs/2];
          if(maxBr<imgRes2[indx+ic1+ir1*xs/2])
            maxBr=imgRes2[indx+ic1+ir1*xs/2];
        }
      }

      if( mask_flashes[indx]==255 && imgRes2[indx]<minBr+10 )
      {
        mask_flashes[indx]=0;
        mask_flashesTmp[indx] = 100;
      }
      if( mask_flashes[indx]==0 && imgRes2[indx]>minBr+10 )
      {
        mask_flashes[indx]=255;
        mask_flashesTmp[indx] = 255;
      }
      if( mask_flashes[indx-1]==0 && imgRes2[indx-1]>minBr+25 )
      {
        mask_flashes[indx]=255;
        mask_flashesTmp[indx] = 255;
      }
      if( mask_flashes[indx+1]==0 && imgRes2[indx+1]>minBr+25 )
      {
        mask_flashes[indx]=255;
        mask_flashesTmp[indx] = 255;
      }
    }
  }

  mask_flashes_Purified = mask_flashes;

  if(sv_f) SAVEMATRIX(mask_flashes_Purified,xs/2,ys/2);


  // solving Poisson's Equation
  // optimization factor for succsesive relaxation
  //double fi=2/(1+(1-(1-2*sin(3.1415926/2/15).^2).^2).^0.5); 
  for( numIter=0; numIter<10; numIter++)
  { 

    for( ir=1; ir<ys/2-1; ir++)
    {
      for( ic=1; ic<xs/2-1; ic++)
      {

//        if( 1 )// use acceleration    
          if( mask_flashes[ir*xs/2+ic]!=0 )
          {
            Img_Iter1[ir*xs/2+ic] = 0.25*( Img_Iter2[(ir-1)*xs/2+ic] + Img_Iter2[ir*xs/2+ic-1] +
              Img_Iter2[(ir+1)*xs/2+ic] + Img_Iter2[ir*xs/2+ic+1] );
            Img_Iter2[ir*xs/2+ic] = fi* Img_Iter1[ir*xs/2+ic] + (1-fi)*Img_Iter2[ir*xs/2+ic];

          }
//
  //        else
    //        if( mask_flashes[ir*xs/2+ic]!=0 )
      //        Img_Iter2[ir*xs/2+ic] = 0.25*( Img_Iter2[(ir-1)*xs/2+ic] + Img_Iter2[ir*xs/2+ic-1] +
        //      Img_Iter2[(ir+1)*xs/2+ic] + Img_Iter2[ir*xs/2+ic+1] );

      }
    }

  }



  // purify mask excuding pixel which are brighter after diffusion
  for( ir=2; ir<ys/2-2; ir++) /* KA: fix from for( ir=2; ir<ys/2-2; ir++), Dr.Phan tested */
  {
    for( ic=2; ic<xs/2-2; ic++)
    {
      int indx = ir*xs/2+ic;
      if( mask_flashes[indx]==255 && Img_Iter2[indx]>imgRes2[indx]+5 && imgRes2[indx]<100 )
      {
        mask_flashes[indx]=0;
        mask_flashesTmp[indx] = 255;
      }
    }
  }

  img_FirstIter = imPnt;
  mask_flashes_ExcBr = mask_flashes;

  if(sv_f) SAVEMATRIX(img_FirstIter,xs,ys);
  if(sv_f) SAVEMATRIX(mask_flashes_ExcBr,xs/2,ys/2);


  // new iteration for difusion
  for( numIter=0; numIter<10; numIter++)
    for( ir=1; ir<ys/2-1; ir++)
      for( ic=1; ic<xs/2-1; ic++)
      {
        if( mask_flashes[ir*xs/2+ic]!=0 )
        {
          Img_Iter1[ir*xs/2+ic] = 0.25*( Img_Iter2[(ir-1)*xs/2+ic] + Img_Iter2[ir*xs/2+ic-1] +
            Img_Iter2[(ir+1)*xs/2+ic] + Img_Iter2[ir*xs/2+ic+1] );
          Img_Iter2[ir*xs/2+ic] = fi* Img_Iter1[ir*xs/2+ic] + (1-fi)*Img_Iter2[ir*xs/2+ic];

        }
      }


      //if(sv_f) SAVEMATRIX(mask_flashes,xs/2,ys/2);



      // copy result to image (whith resizing) 
      for( ir=1; ir<ys/2-1; ir++)
      {
        for( ic=1; ic<xs/2-1; ic++)
        {
          if( mask_flashes[ir*xs/2+ic]!=0 )
          {
            double tmpd = Img_Iter2[ir*xs/2+ic];
            if(tmpd<0) tmpd = 0;
            if(tmpd>255) tmpd = 255;
            //if(imPnt[ir*2*xs+ic*2] > tmpd)
            imPnt[ir*2*xs+ic*2] = (unsigned char)tmpd;
            //if(imPnt[ir*2*xs+ic*2+1] > tmpd)
            imPnt[ir*2*xs+ic*2+1] = (unsigned char)tmpd;
            //if(imPnt[(ir*2+1)*xs+ic*2] > tmpd)
            imPnt[(ir*2+1)*xs+ic*2] = (unsigned char)tmpd;
            //if(imPnt[(ir*2+1)*xs+ic*2+1] > tmpd)
            imPnt[(ir*2+1)*xs+ic*2+1] = (unsigned char)tmpd;

            // mitigate interlasing effect
            if(imPnt[ir*2*xs+ic*2-1] > tmpd+10)
              imPnt[ir*2*xs+ic*2-1] = (unsigned char)((imPnt[ir*2*xs+ic*2-1]+tmpd)/2);
            if(imPnt[ir*2*xs+ic*2+2] > tmpd+10)
              imPnt[ir*2*xs+ic*2+2] = (unsigned char)((imPnt[ir*2*xs+ic*2+2]+tmpd)/2);
            if(imPnt[(ir*2+1)*xs+ic*2-1] > tmpd+10)
              imPnt[(ir*2+1)*xs+ic*2-1] = (unsigned char)((imPnt[(ir*2+1)*xs+ic*2-1]+tmpd)/2);
            if(imPnt[(ir*2+1)*xs+ic*2+2] > tmpd+10)
              imPnt[(ir*2+1)*xs+ic*2+2] = (unsigned char)((imPnt[(ir*2+1)*xs+ic*2+2]+tmpd)/2);

          }

//          if(0)
  //        {
    //        imPnt[ir*2*xs+ic*2] = mask_flashes[ir*xs/2+ic];
      //      imPnt[ir*2*xs+ic*2+1] = mask_flashes[ir*xs/2+ic];
        //    imPnt[(ir*2+1)*xs+ic*2] = mask_flashes[ir*xs/2+ic];
          //  imPnt[(ir*2+1)*xs+ic*2+1] = mask_flashes[ir*xs/2+ic];			
          //}
        }
      }

      img_Final = imPnt;
      mask_flashes_Final = mask_flashes;
      if(sv_f) SAVEMATRIX(img_Final,xs,ys);
      if(sv_f) SAVEMATRIX(mask_flashes_Final,xs/2,ys/2);
      // free all memory for processing
      free(mask_flashes0);
      // return DLLT_OK;
  return 0;
}

void IP_IVIR_PMP_Erode3x3(
              unsigned char* dst,
              const unsigned char* src,
              int xs,
              int ys)
{
  int i,j,m,q;

  // upper boundary
  memset(dst,0,xs);
  // body
  for (i=1;i<ys-1;i++)
  {
    // left
    dst[i*xs] = 0;
    // body
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
    // right
    dst[i*xs+xs-1] = 0;
  }
  // lower boundary
  memset(dst+(ys-1)*xs,0,xs);
}


void IP_IVIR_PMP_Erode3x3Cross(
                 unsigned char* dst,
                 const unsigned char* src,
                 int xs,
                 int ys)
{
  int i,j,m,q;

  // upper boundary
  memset(dst,0,xs);
  // body
  for (i=1;i<ys-1;i++)
  {
    // left
    dst[i*xs] = 0;
    // body
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
    // right
    dst[i*xs+xs-1] = 0;
  }
  // lower boundary
  memset(dst+(ys-1)*xs,0,xs);
}


void IP_IVIR_PMP_Dilate3x3(
               unsigned char* dst,
               const unsigned char* src,
               int xs,
               int ys)
{
  int i,j,m,q;

  // upper boundary
  memset(dst,0,xs);
  // body
  for (i=1;i<ys-1;i++)
  {
    // left
    dst[i*xs] = 0;
    // body
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
    // right
    dst[i*xs+xs-1] = 0;
  }
  // lower boundary
  memset(dst+(ys-1)*xs,0,xs);
}

void IP_IVIR_PMP_Dilate3x3Cross(
                unsigned char* dst,
                const unsigned char* src,
                int xs,
                int ys)
{
  int i,j,m,q;

  // upper boundary
  memset(dst,0,xs);
  // body
  for (i=1;i<ys-1;i++)
  {
    // left
    dst[i*xs] = 0;
    // body
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
    // right
    dst[i*xs+xs-1] = 0;
  }
  // lower boundary
  memset(dst+(ys-1)*xs,0,xs);
}

void IP_IVIR_PMP_FillBorders1B(
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



// remove flashes by harmonic filling
int IP_IVIR_DetectFlashes(
              unsigned char* im,    // IN:  input image
              int xs,                     // IN:  image size
              int ys,
              unsigned char* mask)    // OUT:  mask
{
  int ir, ic, tmpint;
  int iAccum, threshDiff1, threshDiff2, threshDiff12;
//	double fi=1.65;		// optimization factor
  int sv_f = 0;      // save or not results
  int i,j;


  //reserve memory for processing of data
  int histVector[256];                   // histogram array for threshhold calculation
  unsigned char* imPnt = &im[0];

  unsigned char* mask_flashes0;
  unsigned char* mask_flashes;
  unsigned char* Im_eroded_1;
  unsigned char* Im_eroded_2;
  unsigned char* Im_dilated_1;
  unsigned char* Im_dilated_2;
  unsigned char* Im_diff;
  unsigned char* imgRes2;  // for reduced image in 2 times
  double* Img_Iter1; // for iteration  during solution poison equation
  double* Img_Iter2; // for itaration durin solution poison equation
  unsigned char* mask_flashesTmp;

  unsigned char* mask_flashes_Purified;

  // optimized reservation
  mask_flashes0 = (unsigned char*) malloc(xs/2*ys/2*8 + xs/2*ys/2*3*8);
  if (!mask_flashes0)
  {
    return -1;
  }
  mask_flashes = (unsigned char*) (mask_flashes0+xs/2*ys/2);
  Im_eroded_1 = (unsigned char*) (mask_flashes+xs/2*ys/2);
  Im_eroded_2 = (unsigned char*) (Im_eroded_1+xs/2*ys/2);
  Im_dilated_1 = (unsigned char*) (Im_eroded_2+xs/2*ys/2);
  Im_dilated_2 = (unsigned char*) (Im_dilated_1+xs/2*ys/2);
  Im_diff = (unsigned char*) (Im_dilated_2+xs/2*ys/2);
  imgRes2 = (unsigned char*) (Im_diff+xs/2*ys/2 );  // for reduced image in 2 times
  Img_Iter1 = (double*) (imgRes2+xs/2*ys/2 ); // for iteration  during solution poison equation
  Img_Iter2 = (double*) (Img_Iter1+xs/2*ys/2 ); // for itaration durin solution poison equation
  mask_flashesTmp = (unsigned char*) (Img_Iter2+xs/2*ys/2);

  // reduce size 2 times
  for( ir=0; ir<ys/2; ir++)
  {
    for( ic=0; ic<xs/2; ic++)
    {
      imgRes2[ir*xs/2+ic] = (unsigned char)( 0.25*( im[ir*2*xs+ic*2] + im[ir*2*xs+ic*2+1] + 
        im[(ir*2+1)*xs+ic*2] + im[(ir*2+1)*xs+ic*2+1] + 2 ) );
      //imgRes2[ir*xs/2+ic] = __min(
      //	__min( im[ir*2*xs+ic*2], im[ir*2*xs+ic*2+1]),
      //	__min(im[(ir*2+1)*xs+ic*2], im[(ir*2+1)*xs+ic*2+1] ) 
      //	) ;
      Img_Iter2[ir*xs/2+ic] = imgRes2[ir*xs/2+ic];
    }
  }
  //if(sv_f) SAVEMATRIX(im,xs,ys);
  if(sv_f) SAVEMATRIX(imgRes2,xs/2,ys/2);


  // erosion-dilation

  // erosion 5 times
  memcpy(Im_eroded_1, imgRes2, xs/2*ys/2);
  for(ic=0; ic<2; ic++)
  {
    IP_IVIR_PMP_Erode3x3(Im_eroded_2,Im_eroded_1,xs/2,ys/2);
    IP_IVIR_PMP_FillBorders1B(  Im_eroded_2,xs/2,ys/2);     // removing borders gap
    memcpy(Im_eroded_1, Im_eroded_2, xs/2*ys/2);
    IP_IVIR_PMP_Erode3x3Cross(Im_eroded_2,Im_eroded_1,xs/2,ys/2);
    IP_IVIR_PMP_FillBorders1B(  Im_eroded_2,xs/2,ys/2);     // removing borders gap
    memcpy(Im_eroded_1, Im_eroded_2, xs/2*ys/2);
  }
  if(sv_f) SAVEMATRIX(Im_eroded_1,xs/2,ys/2);

  // dilation 5 times
  memcpy(Im_dilated_1, Im_eroded_1, xs/2*ys/2);
  for(ic=0; ic<2; ic++)
  {
    IP_IVIR_PMP_Dilate3x3Cross(Im_dilated_2,Im_dilated_1,xs/2,ys/2);
    memcpy(Im_dilated_1, Im_dilated_2, xs/2*ys/2);
    IP_IVIR_PMP_Dilate3x3(Im_dilated_2,Im_dilated_1,xs/2,ys/2);
    memcpy(Im_dilated_1, Im_dilated_2, xs/2*ys/2);
  }
  if(sv_f) SAVEMATRIX(Im_dilated_1,xs/2,ys/2);


  // calculate opening: difference of image and opening
  // and fill histogramm
  memset(histVector, 0, 256*sizeof(int) );
  for( ir=0; ir<ys/2*xs/2; ir++)
  {
    tmpint = imgRes2[ir]-Im_eroded_2[ir]; 
    if(tmpint<0) tmpint = 0; 
    if(tmpint>255) tmpint = 255; 

    Im_diff[ir] = (unsigned char)tmpint;
    //if()
    histVector[Im_diff[ir]]++; 
  }
  if(sv_f) SAVEMATRIX(Im_diff,xs/2,ys/2);

  // calculate threshhold as 
  //threshDiff1 = prctile(imgDiff(:), 100.*(1.-50/prod(size(imgDiff)) )  );
  //threshDiff2 = prctile(imgDiff(:), 99.5  );
  //threshDiff2 = threshDiff1-10;
  //threshDiff12 = ( threshDiff1+2*threshDiff2)/3;

  iAccum = 0;
  for( ir=255; ir>0; ir--)
  {
    iAccum+=histVector[ir];
    if(iAccum>50) // 50 pixels from top
      break;
  }
  threshDiff1 = ir;

  iAccum = 0;
  for( ir=255; ir>0; ir--)
  {
    iAccum+=histVector[ir];
    if(iAccum>ys/2*xs/2/200) // 99.5%
      break;
  }
  threshDiff2 = ir-10;
  threshDiff12 = ( threshDiff1+2*threshDiff2)/3;

  // use threshold for masking
  for( ir=0; ir<ys/2*xs/2; ir++)
  {
    mask_flashes0[ir] = 0;
    if( Im_diff[ir]>threshDiff12 )
      mask_flashes0[ir] = 255;
  }
  //if(sv_f) ; VCL_BMP_Save(mask_flashes0,xs/2,ys/2,xs/2,"c:\\Im_t_masked_0.bmp",USUAL_BMP);

  // dilate mask
  IP_IVIR_PMP_Dilate3x3Cross(mask_flashes, mask_flashes0,xs/2,ys/2);
  memcpy(mask_flashes0, mask_flashes, xs/2*ys/2);
  IP_IVIR_PMP_Dilate3x3(mask_flashes, mask_flashes0,xs/2,ys/2);
  memcpy(mask_flashes0, mask_flashes, xs/2*ys/2);
  IP_IVIR_PMP_Dilate3x3Cross(mask_flashes, mask_flashes0,xs/2,ys/2);

  if(sv_f) SAVEMATRIX(mask_flashes0,xs/2,ys/2);


  // tune for brightness
  for( ir=0; ir<ys/2*xs/2; ir++)
  {
    mask_flashesTmp[ir] = 0;
    if( imgRes2[ir]>200 )
      mask_flashesTmp[ir] = 255;
  }

  // dilate mask
  IP_IVIR_PMP_Dilate3x3(mask_flashes0, mask_flashesTmp,xs/2,ys/2);
  memcpy(mask_flashesTmp, mask_flashes0, xs/2*ys/2);
  IP_IVIR_PMP_Dilate3x3Cross(mask_flashes0, mask_flashesTmp,xs/2,ys/2);
  memcpy(mask_flashesTmp, mask_flashes0, xs/2*ys/2);
  IP_IVIR_PMP_Dilate3x3(mask_flashes0, mask_flashesTmp,xs/2,ys/2);
  memcpy(mask_flashesTmp, mask_flashes0, xs/2*ys/2);
  IP_IVIR_PMP_Dilate3x3Cross(mask_flashes0, mask_flashesTmp,xs/2,ys/2);
  memcpy(mask_flashesTmp, mask_flashes0, xs/2*ys/2);
  IP_IVIR_PMP_Dilate3x3(mask_flashes0, mask_flashesTmp,xs/2,ys/2);
  memcpy(mask_flashesTmp, mask_flashes0, xs/2*ys/2);
  IP_IVIR_PMP_Dilate3x3Cross(mask_flashes0, mask_flashesTmp,xs/2,ys/2);

  for( ir=0; ir<ys/2*xs/2; ir++)
    mask_flashes[ir] = __min( mask_flashes0[ir], mask_flashes[ir] );

  if(sv_f) SAVEMATRIX(mask_flashes,xs/2,ys/2);


  // correct mask to exclude minimal brightness
  for( ir=3; ir<ys/2-3; ir++)
  {
    for( ic=3; ic<xs/2-3; ic++)
    {
      int isVertMaskBr;
      int isHorzMaskBr;
      int isBr;
      double minBr;
      double maxBr;
      int ir1;
      int ic1;
      int indx = ir*xs/2+ic;
      mask_flashesTmp[indx] = 0;

      isVertMaskBr = mask_flashes[indx]+mask_flashes[indx-xs/2]+mask_flashes[indx+xs/2];
      isHorzMaskBr = mask_flashes[indx]+mask_flashes[indx-1]+mask_flashes[indx+1];
      isBr = 0;
      isVertMaskBr /= 255;
      isHorzMaskBr /= 255;
      if(isVertMaskBr==0 || isVertMaskBr==3)
        continue;
      if(isHorzMaskBr==0 || isHorzMaskBr==3)
        continue;

      minBr = 10000;
      maxBr = 0;
      for(ir1=-2;ir1<=2;ir1++)
      {
        for(ic1=-2;ic1<=2;ic1++)
        {
          if(minBr>imgRes2[indx+ic1+ir1*xs/2])
            minBr=imgRes2[indx+ic1+ir1*xs/2];
          if(maxBr<imgRes2[indx+ic1+ir1*xs/2])
            maxBr=imgRes2[indx+ic1+ir1*xs/2];
        }
      }

      if( mask_flashes[indx]==255 && imgRes2[indx]<minBr+10 )
      {
        mask_flashes[indx]=0;
        mask_flashesTmp[indx] = 100;
      }
      if( mask_flashes[indx]==0 && imgRes2[indx]>minBr+10 )
      {
        mask_flashes[indx]=255;
        mask_flashesTmp[indx] = 255;
      }
      if( mask_flashes[indx-1]==0 && imgRes2[indx-1]>minBr+25 )
      {
        mask_flashes[indx]=255;
        mask_flashesTmp[indx] = 255;
      }
      if( mask_flashes[indx+1]==0 && imgRes2[indx+1]>minBr+25 )
      {
        mask_flashes[indx]=255;
        mask_flashesTmp[indx] = 255;
      }
    }
  }

  mask_flashes_Purified = mask_flashes;

  if(sv_f) SAVEMATRIX(mask_flashes_Purified,xs/2,ys/2);

  imPnt = mask;

  // copy result to image (whith resizing) 
  for( ir=1; ir<ys/2-1; ir++)
  {
    for( ic=1; ic<xs/2-1; ic++)
    {
      if( mask_flashes[ir*xs/2+ic]!=0 )
      {
        double tmpd = mask_flashes_Purified[ir*xs/2+ic];
        if(tmpd<0) tmpd = 0;
        if(tmpd>255) tmpd = 255;
        //if(imPnt[ir*2*xs+ic*2] > tmpd)
        imPnt[ir*2*xs+ic*2] = (unsigned char)tmpd;
        //if(imPnt[ir*2*xs+ic*2+1] > tmpd)
        imPnt[ir*2*xs+ic*2+1] = (unsigned char)tmpd;
        //if(imPnt[(ir*2+1)*xs+ic*2] > tmpd)
        imPnt[(ir*2+1)*xs+ic*2] = (unsigned char)tmpd;
        //if(imPnt[(ir*2+1)*xs+ic*2+1] > tmpd)
        imPnt[(ir*2+1)*xs+ic*2+1] = (unsigned char)tmpd;

        // mitigate interlasing effect
        if(imPnt[ir*2*xs+ic*2-1] > tmpd+10)
          imPnt[ir*2*xs+ic*2-1] = (unsigned char)((imPnt[ir*2*xs+ic*2-1]+tmpd)/2);
        if(imPnt[ir*2*xs+ic*2+2] > tmpd+10)
          imPnt[ir*2*xs+ic*2+2] = (unsigned char)((imPnt[ir*2*xs+ic*2+2]+tmpd)/2);
        if(imPnt[(ir*2+1)*xs+ic*2-1] > tmpd+10)
          imPnt[(ir*2+1)*xs+ic*2-1] = (unsigned char)((imPnt[(ir*2+1)*xs+ic*2-1]+tmpd)/2);
        if(imPnt[(ir*2+1)*xs+ic*2+2] > tmpd+10)
          imPnt[(ir*2+1)*xs+ic*2+2] = (unsigned char)((imPnt[(ir*2+1)*xs+ic*2+2]+tmpd)/2);

      }
    }
  }

  if(sv_f) SAVEMATRIX(mask,xs,ys);

  free(mask_flashes0);

  return 0;
}
