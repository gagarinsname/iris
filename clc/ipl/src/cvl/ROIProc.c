/**********************************************************************    
  ROIProc.c
 
   Project
      Computer Vision Library (CVL)

   Description
       Implementation of functions, used to process with ROI.

   Author
      Alexander Gneushev

 **********************************************************************/
#include "ipl.h"

#define max(a,b) (a>b?a:b)
#define min(a,b) (a>b?b:a)

// ExpandRectROI function
// expandes rect ROI on 2*HalfWidthSizeExpand and 2*HalfHeightSizeExpand values
void CVL_ExpandRectROI(
  int*                pDestROIStartX,        // Dest ROI left
  int*                pDestROIStartY,        // Dest ROI top
  unsigned int*       pDestROIWidth,         // Dest ROI width  
  unsigned int*       pDestROIHeight,        // Dest ROI height
  int                 BoundROIStartX,        // Bounding ROI left
  int                 BoundROIStartY,        // Bounding ROI top
  unsigned int        BoundROIWidth,         // Bounding ROI width 
  unsigned int        BoundROIHeight,        // Bounding ROI height
  const unsigned int  HalfWidthSizeExpand,   // Half of width expand value
  const unsigned int  HalfHeightSizeExpand)  // Half of height expand value
{
  *pDestROIStartX -= (int) HalfWidthSizeExpand;
  *pDestROIWidth  += 2*HalfWidthSizeExpand;
  *pDestROIStartY -= (int) HalfHeightSizeExpand;
  *pDestROIHeight += 2*HalfHeightSizeExpand;

  if (*pDestROIStartX < BoundROIStartX + (int)BoundROIWidth &&
      *pDestROIStartY < BoundROIStartY + (int)BoundROIHeight &&
      *pDestROIStartX + (int)*pDestROIWidth > BoundROIStartX &&
      *pDestROIStartY + (int)*pDestROIHeight > BoundROIStartY)
   {
     if (*pDestROIStartX < BoundROIStartX) 
     {
         *pDestROIWidth -= BoundROIStartX - *pDestROIStartX;
         *pDestROIStartX = BoundROIStartX;
     }
     if (*pDestROIStartX + (int)*pDestROIWidth > BoundROIStartX + (int)BoundROIWidth)
            *pDestROIWidth = BoundROIStartX + (int)BoundROIWidth - *pDestROIStartX;
     if (*pDestROIStartY < BoundROIStartY) 
     {
         *pDestROIHeight -= BoundROIStartY - *pDestROIStartY;
         *pDestROIStartY = BoundROIStartY;
     }
     if (*pDestROIStartY + (int)*pDestROIHeight > BoundROIStartY + (int)BoundROIHeight)
            *pDestROIHeight = BoundROIStartY + (int)BoundROIHeight - *pDestROIStartY;
  }
  else
  {
   *pDestROIStartX = 0;
   *pDestROIStartY = 0;
   *pDestROIWidth  = 0;
   *pDestROIHeight = 0;
  }

}
