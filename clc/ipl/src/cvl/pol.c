/**********************************************************************    
  pol.c
 
   Project
      Computer Vision Library (CVL)

   Author
      Alexander Gneushev

 **********************************************************************/
#include <math.h>
#include "bpl.h"
#include "ipl.h"

void IPL_POLAR_TransformBilinearHorFillByLast(
  unsigned char*       Dest, 
  const unsigned char* Source,
  unsigned int         DestWidth,
  unsigned int         DestHeight,
  unsigned int         DestWidthItems,
  unsigned int         SourceWidth,
  unsigned int         SourceHeight,
  unsigned int         SourceWidthItems,
  int                  SrcCx, 
  int                  SrcCy,
  unsigned int         SrcStartRadius,
  unsigned int         SrcFinishRadius,
  float                SrcStartAngleRad,
  float                SrcFinishAngleRad,
  const unsigned char  FirstFillValue,
  int*                 HorizontalBuffer,    // SizeItems=DestWidth
  float*               VerticalBuffer)      // SizeItems=DestHeight
{
  unsigned char*       pDest;
  const unsigned char* pSource;
  unsigned int AddToDest = DestWidthItems - DestWidth;
  unsigned int  i,j;
  unsigned char Value;
  int           x,y,xx,yy,kx,ky;
  int           Cosa,Sina;
  const int          Cx = 256*SrcCx;
  const int          Cy = 256*SrcCy;
  const float        AngleRange  = SrcFinishAngleRad-SrcStartAngleRad;
  const unsigned int RadiusRange = SrcFinishRadius-SrcStartRadius;

  float* yj = VerticalBuffer;
  int*   xi = HorizontalBuffer;

  if( Dest &&  (DestWidth>1) && (DestHeight>1) && (DestWidth<= DestWidthItems)&&
      Source && (SourceWidth>1) && (SourceHeight>1) && (SourceWidth<= SourceWidthItems)&&
      ((SrcCx+(int)SrcFinishRadius) > 0)&&((SrcCy+(int)SrcFinishRadius)>0)&&
      ((SrcCx-(int)SrcFinishRadius) < (int)SourceWidth)&&
      ((SrcCy-(int)SrcFinishRadius) < (int)SourceHeight)&&
      (SrcStartRadius < SrcFinishRadius)&&(fabs(AngleRange)>0.0f))
  {	
    pDest = Dest;
    
    for (j = 0; j<DestHeight; j++)
      yj[j] = SrcStartAngleRad+(j*AngleRange)/DestHeight;
  
    for (i=0;i<DestWidth;i++)
        xi[i] = (int)(256*SrcStartRadius+(256*i*RadiusRange)/DestWidth);

    for (j = 0; j<DestHeight; j++)
    {     
//      Cosa   = (int)(cos(yj[j])*0x1000);     
      Cosa = (MIA_FIX20_cos((int)(yj[j]*(1<<20))))/(1<<8);
//      Sina   = (int)(sin(yj[j])*0x1000);
      Sina = (MIA_FIX20_sin((int)(yj[j]*(1<<20))))/(1<<8);
      Value  = FirstFillValue; 

      for(i = 0; i < DestWidth; i++)
      {       
        xx = Cx + ((xi[i]*Cosa)>>12);
        yy = Cy - ((xi[i]*Sina)>>12);
        
        x  = xx>>8;
        y  = yy>>8;
        ky = yy&255;
        kx = xx&255;
               
        if ((x>=0)&&(x < (int)SourceWidth-1)&&(y >= 0)&&(y < (int)SourceHeight-1))
        {
          pSource=Source+x+y*SourceWidthItems;
          Value = (unsigned char)((((pSource[SourceWidthItems+1]*kx + (pSource[SourceWidthItems])*(256-kx)))*ky + ((pSource[1]*kx + (pSource[0])*(256-kx)))*(256-ky))>>16);
        }
       
       *pDest++ = Value;
      }

      pDest += AddToDest;
    }    
  }  	
}

void IPL_POLAR_MakePolarRubberInt(
  uint8* polar,
  int32 w,
  int32 numLines,
  const uint8* im,
  int32 xs,
  int32 ys,
  int32 xC, 
  int32 yC, 
  int32 rC, 
  int32 xCI, 
  int32 yCI, 
  int32 rCI)
{
  int32 i,angleInt,idx,Cosa12,Sina12;
  int32 Pi_FIX20 = 3294199, angle_FIX20, x1_FIX12, x2_FIX12, y1_FIX12, y2_FIX12, x, y;
  //double x,y,angle,x1,x2,y1,y2;
  //double piMultiplier = 2.*3.1415926/numLines;
  int32 R = w, R_FIX12 = w<<12;
  uint8 imBr=0;

  for(angleInt=0;angleInt<numLines;angleInt++)
  {
    angle_FIX20 = Pi_FIX20*angleInt*2+Pi_FIX20/2+1;
    Cosa12 = MIA_FIX20_cos(angle_FIX20)>>8;
    Sina12 = MIA_FIX20_sin(angle_FIX20)>>8;

    x1_FIX12 = Cosa12*rC  + (xC<<12);
    x2_FIX12 = Cosa12*rCI + (xCI<<12);
    y1_FIX12 =-Sina12*rC/100*94  + (yC<<12);
    y2_FIX12 =-Sina12*rCI/100*95 + (yCI<<12);

    for( i=0;i<R;i++ )
    {
      x = ( (x1_FIX12*(R-i) + x2_FIX12*i + R_FIX12/2) /R )>>12;
      y = ( (y1_FIX12*(R-i) + y2_FIX12*i + R_FIX12/2) /R )>>12;
      if( y<1 || y>=ys-1 || x<1 || x>=xs-1 )
      {
        polar[ angleInt*w+i ] = imBr;
        continue;
      }
      idx = (x+y*xs);
      polar[ angleInt*w+i ] = im[idx];
      imBr = im[idx];
    }
  }
}

void IPL_POLAR_MakePolarRubber(
  uint8* polar,
  int32 w,
  int32 numLines,
  const uint8* im,
  int32 xs,
  int32 ys,
  int32 xC, 
  int32 yC, 
  int32 rC, 
  int32 xCI, 
  int32 yCI, 
  int32 rCI)
{
  int32 i,angleInt,idx;
  double x,y,angle, x1,x2,y1,y2, cosv, sinv;
  const double piMultiplier = 2.*3.141592653589/numLines;
  const double anglShift = (3.141592653589/2);
  int32 R = w;
  uint8 imBr=0;

  for(angleInt=0;angleInt<numLines;angleInt++)
  {
    angle = piMultiplier*angleInt+anglShift;
    cosv = cos(angle);
    sinv = sin(angle);
    x1=cosv*rC+xC;
    x2=cosv*rCI+xCI;
    y1=-sinv*rC*.94+yC;
    y2=-sinv*rCI*.95+yCI;
    for( i=0;i<R;i++ )
    {
      x = (int)( (x1*(R-i) + x2*i)/R + 0.5);
      y = (int)( (y1*(R-i) + y2*i)/R + 0.5);
      if( y<1 || y>=ys-1 || x<1 || x>=xs-1 )
      {
        polar[ angleInt*w+i ] = imBr;
        continue;
      }
      idx = (int)(x+y*xs);
      polar[ angleInt*w+i ] = im[idx];
      imBr = im[idx];
    }
  }
}

// draw polar lines on original image
void IPL_POLAR_DrawPolarRubberLines(
  uint8* im,           // IN:OUT
  int32 xs,
  int32 ys,
  int32 xC, 
  int32 yC, 
  int32 rC, 
  int32 xCI, 
  int32 yCI, 
  int32 rCI)
{
  int32 i,angleInt,idx;
  double x,y,angle, x1,x2,y1,y2, cosv, sinv;
  const int32 numLines = 256;  // number of angular 
  const int32 w = 32;  // number of angular 
  const double piMultiplier = 2.*3.141592653589/numLines;
  const double anglShift = (3.141592653589/2);
  int32 R = w;

    for(angleInt=0;angleInt<numLines;angleInt++)
    {
      angle = piMultiplier*angleInt+anglShift;
      cosv = cos(angle);
      sinv = sin(angle);
      x1=cosv*rC+xC;
      x2=cosv*rCI+xCI;
      y1=-sinv*rC*.94+yC;
      y2=-sinv*rCI*.95+yCI;
      for( i=0;i<R;i+=4 )
      {
        x = (int)( (x1*(R-i) + x2*i)/R + 0.5);
        y = (int)( (y1*(R-i) + y2*i)/R + 0.5);
        if( y<1 || y>=ys-1 || x<1 || x>=xs-1 )
        {
          continue;
        }
        idx = (int)(x+y*xs);
        im[idx]=255;
      }
    }
}

void IPL_POLAR_MakePolarRubberDeform(
  uint8* polar,
  int32 w,
  int32 numLines,
  const uint8* im,
  int32 xs,
  int32 ys,
  int32 xC, 
  int32 yC, 
  int32 rC, 
  int32 xCI, 
  int32 yCI, 
  int32 rCI,
  int32 DeformMethod, // deformation method; 0-usual rubber; 1-power mode; 2-Ivan deformation 
  double dscaler)     // parameter of deformation for method 1 and 2
{
  int32 i,angleInt,idx;
  double x,y,angle, x1,x2,y1,y2, cosv, sinv, dX, dY, i_R, i_R_pow, newV;
  const double piMultiplier = 2.*3.141592653589/numLines;
  const double anglShift = (3.141592653589/2);
  int32 R = w;
  uint8 imBr=0;
  //for Ivan deformation
  double r_p = 1.*rC/rCI;
  double alfa = .5, ro_dd, ro_d, ro, SS, U, V, x_d, y_d;
  double tetta1 = 0.95, tetta2 = 1;
  const int32 DRAW_LINE = 0;   // draw epipolar lines on image
                               // only for visualization - vialite template information

  
  if( alfa==0 && DeformMethod==2)
    alfa = r_p;
  else
    alfa = dscaler;

  //alfa = r_p;

  if(DeformMethod<0 || DeformMethod>2)
    DeformMethod=0;
  for(angleInt=0;angleInt<numLines;angleInt++)
  {
    angle = piMultiplier*angleInt+anglShift;
    cosv = cos(angle);
    sinv = sin(angle);
    x1=cosv*rC+xC;
    x2=cosv*rCI+xCI;
    y1=-sinv*rC*.94+yC;
    y2=-sinv*rCI*.95+yCI;
    for( i=0;i<R;i++ )
    {
      if(DeformMethod==0)
      {
        x = (int)( (x1*(R-i) + x2*i)/R + 0.5);
        y = (int)( (y1*(R-i) + y2*i)/R + 0.5);
        if( y<1 || y>=ys-1 || x<1 || x>=xs-1 )
        {
          polar[ angleInt*w+i ] = imBr;
          continue;
        }
        idx = (int)(x+y*xs);
        if(DRAW_LINE && ((i%4)==0 || i==R-1) ) // draw epipolar points on original image
        {
          uint8* pimDraw = (uint8*)(im+idx);
          *pimDraw = 255;
          if(angleInt%2) *pimDraw = 0;
        }        
        polar[ angleInt*w+i ] = im[idx];
        imBr = im[idx];      
      }
      if(DeformMethod==1)  // power mode
      {
        //iV=((0:(R-1))./(R-1)).^optionsTMP.tmplgen_polarscaler*(R-1);
        //x = ( (x1*(R-i) + x2*i)/R + 0.5);
        //y = ( (y1*(R-i) + y2*i)/R + 0.5);
        //x = x1*(1.-i/R) + x2*i/R ;
        //y = y1*(1.-i/R) + y2*i/R ;
        i_R = 1.*i/R;
        i_R_pow = pow(i_R, dscaler);
        x = x1*(1.-i_R_pow)+x2*i_R_pow;
        y = y1*(1.-i_R_pow)+y2*i_R_pow;
        if( y<1 || y>=ys-1 || x<1 || x>=xs-1 )
        {
          polar[ angleInt*w+i ] = imBr;
          continue;
        }

        // bilinear interpolation
			  //s = xi->data[i]-floor(xi->data[i]);
			  //t = yi->data[i]-floor(yi->data[i]);
        //zi->data[i] = (z->data[ndx]*(1-t)+z->data[ndx+ncols]*t)*(1-s)+(z->data[ndx+1]*(1-t)+z->data[ndx+ncols+1]*t)*s;
        dX = x-(int)(x);
        dY = y-(int)(y);
        idx = (int)(x)+(int)(y)*xs;
        if(DRAW_LINE && ((i%4)==0 || i==R-1) ) // draw epipolar points on original image
        {
          uint8* pimDraw = (uint8*)(im+idx);
          *pimDraw = 255;
          if(angleInt%2) *pimDraw = 0;
        }        
        newV = (im[idx]*(1-dY)+im[idx+xs]*dY)*(1-dX) + (im[idx+1]*(1-dY)+im[idx+xs+1]*dY)*dX;
        polar[ angleInt*w+i ] = (uint8)((int)(newV+.5));
        imBr = polar[ angleInt*w+i ];
      }
      if(DeformMethod==2)  // new Ivan method
      {
        ro_dd = pow(1.*i/R,1.5);
        ro_d = alfa+(1-alfa)*ro_dd;
    
        SS = (alfa-r_p)*r_p/(1-r_p);
        //ro = (ro_d+SS+((ro_d+SS).^2-4*SS).^.5 )/2;
        ro = (ro_d+SS+sqrt((ro_d+SS)*(ro_d+SS)-4*SS) )/2;

        U = ro*cosv;
        V = ro*sinv;
        
        x_d = U*rCI;
        y_d = V*rCI;

        x = xC-(rC*cosv+(x_d/cosv-rC)*(cosv+(xC-xCI)/(rCI-rC) ));
        y = yC-(rC*sinv+(y_d/sinv-rC)*(sinv*(rCI*tetta1-rC)/(rCI*tetta2-rC)+(yC-yCI)/(rCI*tetta2-rC) ));

        if( y<1 || y>=ys-1 || x<1 || x>=xs-1 )
        {
          polar[ angleInt*w+i ] = imBr;
          continue;
        }
        // bilinear interpolation
			  //s = xi->data[i]-floor(xi->data[i]);
			  //t = yi->data[i]-floor(yi->data[i]);
        //zi->data[i] = (z->data[ndx]*(1-t)+z->data[ndx+ncols]*t)*(1-s)+(z->data[ndx+1]*(1-t)+z->data[ndx+ncols+1]*t)*s;
        dX = x-(int)(x);
        dY = y-(int)(y);
        idx = (int)(x)+(int)(y)*xs;
        if(DRAW_LINE && ((i%4)==0 || i==R-1) ) // draw epipolar points on original image
        {
          uint8* pimDraw = (uint8*)(im+idx);
          *pimDraw = 255;
          if(angleInt%2) *pimDraw = 0;
        }
        newV = (im[idx]*(1-dY)+im[idx+xs]*dY)*(1-dX) + (im[idx+1]*(1-dY)+im[idx+xs+1]*dY)*dX;
        polar[ angleInt*w+i ] = (uint8)((int)(newV+.5));
        imBr = polar[ angleInt*w+i ]; 
      }
    }
  }
}

void IPL_POLAR_MakePolarRubberTranspose(
  uint8* polar,
  int32 numLines,
  int32 w,
  const uint8* im,
  int32 xs,
  int32 ys,
  int32 xC, 
  int32 yC, 
  int32 rC, 
  int32 xCI, 
  int32 yCI, 
  int32 rCI)
{
  int32 i,angleInt,idx;
  double x,y,angle, x1,x2,y1,y2;
  double piMultiplaer = 2.*3.1415926/numLines;
  int32 R = w;
  uint8 imBr=0;

  for(angleInt=0;angleInt<numLines;angleInt++)
  {
    angle = piMultiplaer*angleInt+3.1415928/2;
    x1=cos(angle)*rC+xC;
    x2=cos(angle)*rCI+xCI;
    y1=-sin(angle)*rC*.94+yC;
    y2=-sin(angle)*rCI*.95+yCI;
    for( i=0;i<R;i++ )
    {
      x = (int)( (x1*(R-i) + x2*i)/R + 0.5);
      y = (int)( (y1*(R-i) + y2*i)/R + 0.5);
      if( y<1 || y>=ys-1 || x<1 || x>=xs-1 )
      {
        polar[ angleInt+i*numLines ] = imBr;
        continue;
      }
      idx = (int)(x+y*xs);
      polar[ angleInt+i*numLines ] = im[idx];
      imBr = im[idx];
    }
  }
}
