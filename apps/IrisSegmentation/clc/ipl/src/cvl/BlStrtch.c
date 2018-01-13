/**********************************************************************    
  BlStrtch.c
 
   Project
      Computer Vision Library (CVL)

   Author
      SVY,Alexander Gneushev

 **********************************************************************/
#include <string.h>
#include "ipl.h"

void CVL_BilinearStretchInROIFill(
  unsigned char*       Dest,                // Destination
  const unsigned char* Source,              // Source
  unsigned int         DestWidth,           // Width of Dest image.
  unsigned int         DestHeight,          // Height of Dest image.
  unsigned int         DestWidthBytes,      // Dest image stride
  unsigned int         SourceWidth,         // Width of Source image.
  unsigned int         SourceHeight,        // Height of Source image.
  unsigned int         SourceWidthBytes,    // Source image stride
  int                  SourceROIStartX,     // ROI on SOURCE image. On DEST image it is 0.
  int                  SourceROIStartY,     // ROI on SOURCE image. On DEST image it is 0.
  unsigned int         SourceROIWidth,      // ROI on SOURCE image. On DEST image it is DestWidth.
  unsigned int         SourceROIHeight,     // ROI on SOURCE image. On DEST image it is DestHeight.
  unsigned int*        HorizontalBuffer,    // SizeItems=DestWidth
  unsigned int*        VerticalBuffer)      // SizeItems=DestHeight
{
  unsigned int* xj;
  unsigned int* yi;
  unsigned int  i,j;
  unsigned int  kx,ky;
  unsigned int  x,y;
  const unsigned char* SourceStart;
  const unsigned char* SourcePtr;
  const unsigned char* SourceLinePtr;
  unsigned char*       DestPtr;
  unsigned int         AddToDest;
  const unsigned int   SourceWidthItems=SourceWidthBytes/sizeof(Source[0]);
  const unsigned int   DestWidthItems  =DestWidthBytes/sizeof(Dest[0]);
  const int            SourceROIFinishX=SourceROIStartX+(int)SourceROIWidth;
  const int            SourceROIFinishY=SourceROIStartY+(int)SourceROIHeight;
  int                  SourceShiftX  = -SourceROIStartX;
  int                  SourceShiftY  = -SourceROIStartY;
  unsigned int         DestROIStartX = 0;
  unsigned int         DestROIStartY = 0;
  unsigned int         DestROIWidth  = DestWidth;
  unsigned int         DestROIHeight = DestHeight;
  unsigned int         DestFinishEdgeWidth;
  unsigned int         DestFinishEdgeHeight;

  if ((SourceWidth>1)&&(SourceHeight>1)&&(SourceROIWidth>1)&&(SourceROIHeight>1)&&
      (SourceWidth<=SourceWidthItems)&&
      (SourceROIStartY<(int)(SourceHeight-1))&&(SourceROIFinishY>1)&&
      (SourceROIStartX<(int)(SourceWidth-1))&&(SourceROIFinishX>1)&&
      (DestWidth)&&(DestHeight)&&(DestWidthItems))
  {
    if (SourceROIStartX < 0)
    {
      SourceROIStartX  = 0;
      DestROIStartX    = (SourceShiftX*DestWidth)/(SourceROIWidth-1);
      DestROIWidth    -= DestROIStartX;
    }
    if (SourceROIFinishX > (int)SourceWidth)
      DestROIWidth = ((SourceShiftX+SourceWidth-1)*DestWidth)/(SourceROIWidth-1) - DestROIStartX;

    if (SourceROIStartY < 0)
    {
      SourceROIStartY  = 0;
      DestROIStartY    = (SourceShiftY*DestHeight)/(SourceROIHeight-1);
      DestROIHeight   -= DestROIStartY;
    }
    if (SourceROIFinishY > (int)SourceHeight)
       DestROIHeight = ((SourceShiftY+SourceHeight-1)*DestHeight)/(SourceROIHeight-1) - DestROIStartY;

    DestFinishEdgeWidth  = DestWidth -(DestROIStartX+DestROIWidth);
    DestFinishEdgeHeight = DestHeight-(DestROIStartY+DestROIHeight);

    if (DestROIWidth&& DestROIHeight&&
       (DestFinishEdgeWidth < DestWidth)&&(DestFinishEdgeHeight < DestHeight))
    {
      xj=HorizontalBuffer;
      yi=VerticalBuffer;
      for (j=0;j<DestROIWidth;j++)
        xj[j] = 256*j*(SourceROIWidth-1)/(DestWidth);
      xj[DestROIWidth-1]&=~255;
      for (i=0;i<DestROIHeight;i++)
        yi[i] = 256*i*(SourceROIHeight-1)/(DestHeight);
      yi[DestROIHeight-1]&=~255;

      SourceStart = Source+SourceROIStartX+SourceROIStartY*SourceWidthItems;
      DestPtr     = Dest+DestROIStartX+DestROIStartY*DestWidthItems;
      AddToDest   = DestWidthItems-DestROIWidth;

      for (i=0;i<DestROIHeight;i++)
      {
        y = yi[i];
        ky = y&255;

        SourceLinePtr=SourceStart+(y>>8)*SourceWidthItems;
        for (j=0;j<DestROIWidth;j++)
          {
            x = xj[j];
            kx = x&255;

            SourcePtr = SourceLinePtr + (x>>8);
            *DestPtr++ = (unsigned char) ((((SourcePtr[SourceWidthItems+1]*kx + (SourcePtr[SourceWidthItems])*(256-kx)))*ky + ((SourcePtr[1]*kx + (SourcePtr[0])*(256-kx)))*(256-ky))>>16);
          };
        DestPtr += AddToDest;
      };


      if (DestROIStartX || DestFinishEdgeWidth)
      {
        i         = DestROIHeight;
        DestPtr   = Dest+DestROIStartY*DestWidthItems;
        AddToDest = DestWidthItems-DestWidth;
        do
        {
          j = DestROIStartX;
          if (j)
          {
           // Value = DestPtr[DestROIStartX];
            do
            {
              *DestPtr++ = 0;//Value;
            }while(--j);
          }

          DestPtr += DestROIWidth;          
          j = DestFinishEdgeWidth;
          if (j)
          {
           // Value = DestPtr[-1];
            do
            {
              *DestPtr++ = 0;//Value;
            }while(--j);
          }

          DestPtr += AddToDest;
        }while(--i);
      }
      if (DestROIStartY)
      {
         i           = DestROIStartY;
         DestPtr     = Dest;
         //DestLinePtr = Dest+DestROIStartY*DestWidthItems;
         do
         {
          //memcpy(DestPtr,DestLinePtr,DestWidthItems);
           memset(DestPtr,0,DestWidthItems);
          DestPtr += DestWidthItems;
         }while(--i);
      };

      if (DestFinishEdgeHeight)
      {
        i             = DestFinishEdgeHeight;
        DestPtr       = Dest+DestWidthItems*(DestHeight-DestFinishEdgeHeight);
        //DestLinePtr   = DestPtr-DestWidthItems;
        do
        {
          //memcpy(DestPtr,DestLinePtr,DestWidthItems);
          memset(DestPtr,0,DestWidthItems);
          DestPtr += DestWidthItems;
        }while(--i);
      };
     
    };
  };
}
