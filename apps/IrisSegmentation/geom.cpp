#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "ipl.h"
#include "tables.h"

// cosine with fix20 integer arithmetic
/*int MIA_FIX20_cos(
  int x)
{
  int neg,xn,tMin;

    // cos is an even function
    if (x<0)
      x = -x;
    // cos is cyclic
    if (x>=(int)(((float)PI)*(1<<21)))
      x %= (int)(((float)PI)*(1<<21));
    // cos(x+pi) = -cos(x)
    if (x>=(int)(((float)PI)*(1<<20)))
    {
      x -= (int)(((float)PI)*(1<<20));
      neg = -1;
    }
    else
      neg = 1;
    // cos(x+pi/2) = -cos(pi/2-x)
    if (x>=(int)(((float)PI)*(1<<19)))
    {
      neg = -neg;
      x = (int)(((float)PI)*(1<<20))-x;
    }
    // table
    xn = x>>10;
    tMin = table_cos[xn];
    tMin += ((table_cos[xn+1]-tMin)*(x-(xn<<10)))>>10;
    return (neg>0)?tMin:(-tMin);
}
*/


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
      Cosa   = (int)(cos(yj[j])*0x1000);     
//      Cosa = (MIA_FIX20_cos((int)(yj[j]*(1<<20))))/(1<<8);
      Sina   = (int)(sin(yj[j])*0x1000);
//      Sina = (MIA_FIX20_sin((int)(yj[j]*(1<<20))))/(1<<8);
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

MIA_RESULT_CODE IPL_GEOM_GetPolarImage0(
  const unsigned char *BufImg,// IN:  eye image
  int     SizeX,        // IN: image size
  int     SizeY,        //
  int     SizeStride,        //
  unsigned char  *bufpol,   // IN-OUT: polar data block 
  int            SizePolarA,   // IN: number of lines in Polar template
  int            SizePolarR,   // IN: number of columns in Polar template
  int            SizePolarStride, // IN: number of columns in Polar template
  int            x,            //center coordinates
  int            y,
  int     BegRad,       // OUT:begining radius polar data block 
  float   BegAng,       // IN:  rotation angle for frame
  void*   MemoryBuff,   // IN:  tmp memory buffer
  int*    MemoryBuffSize)
{
  int*   HorizontalBuffer = NULL;
  float* VerticalBuffer   = NULL;
  int iNeedMemSize;

    if (!MemoryBuffSize)
      return ERR_GEN_NULLPOINTER;
    iNeedMemSize = SizePolarR*sizeof(HorizontalBuffer[0]) + SizePolarA*sizeof(VerticalBuffer[0]);
    if (MemoryBuff == NULL)
    {
      *MemoryBuffSize = iNeedMemSize;
      return ERR_OK;
    }
    if (*MemoryBuffSize < iNeedMemSize)
    {
      *MemoryBuffSize = iNeedMemSize;
      return ERR_GEN_NOMEMORY;
    }
    HorizontalBuffer = (int*)MemoryBuff;
    VerticalBuffer  = (float*)((unsigned char*)HorizontalBuffer + SizePolarR*sizeof(HorizontalBuffer[0]));
    if (!HorizontalBuffer && !VerticalBuffer)
      return ERR_GEN_INSUFFICIENT_BUFFER; 
    IPL_POLAR_TransformBilinearHorFillByLast( 
        bufpol, 
        BufImg,
        SizePolarR,
        SizePolarA,
        SizePolarStride,
        SizeX,
        SizeY,
        SizeStride,
        x,
        y,
        BegRad,
        SizePolarR,
        0,
        PI*2,
        0,
        HorizontalBuffer,
        VerticalBuffer); 
    return ERR_OK;
}


