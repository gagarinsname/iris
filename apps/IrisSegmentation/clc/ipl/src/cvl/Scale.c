/*************************************************************
  Scale.h

  Description: Declarations of functions, used to scale images.

  Revision history:

  Author   Reason

  SVY      initial C version
  ANG      revised and enhanced
 *************************************************************/
#include <string.h>
#include "unitypes.h"

// This routine used to scale down Source image in ROI.
// Each point of Dest calculated as average of corresponding 4 points of Source.
// So each Dest(X,Y)=(Source(2*X,2*Y)+Source(2*X+1,2*Y)+Source(2*X,2*Y+1)+Source(2*X+1,2*Y+1))/4.
void CVL_DownScaleImage2TimesByCalcAverageInRect(
  unsigned char*       Dest,              //Destination
  const unsigned char* Source,            //Source
  unsigned int         DestROIStartX,
  unsigned int         DestROIStartY,
  unsigned int         SourceROIStartX,   // Source ROI
  unsigned int         SourceROIStartY,   //
  unsigned int         SourceROIWidth,    //
  unsigned int         SourceROIHeight,   //
  unsigned int         DestWidthBytes,    //Row stride of destination image
  unsigned int         SourceWidthBytes)  //Row stride of source image
{
  unsigned char*       pDest;
  const unsigned char* pSource;
  unsigned int         VrtCounter;
  unsigned int         HrzCounter;
  const unsigned int   SourceWidthItems=SourceWidthBytes/sizeof(Source[0]);
  const unsigned int   DestWidthItems  =DestWidthBytes/sizeof(Dest[0]);
  unsigned int         DestROIWidth   = (SourceROIWidth+1)/2;
  unsigned int         DestROIHeight  = SourceROIHeight/2;
  unsigned int         AddToSource;
  unsigned int         AddToDest;
  /*FASTUINT16*/unsigned short           Value;

  if (SourceROIWidth&&SourceROIHeight&&DestROIWidth&&DestROIHeight&&
      ((SourceROIStartX + SourceROIWidth)<=SourceWidthItems)&&
      ((DestROIStartX + DestROIWidth)<=DestWidthItems))
    {
      AddToSource = 2*SourceWidthItems - SourceROIWidth;
      AddToDest   = DestWidthItems - DestROIWidth;
      pSource     = Source + SourceROIStartX + SourceWidthItems*SourceROIStartY;
      pDest       = Dest+DestROIStartX + DestWidthItems*DestROIStartY;

      VrtCounter  = DestROIHeight;
      if (SourceROIStartX&1) //SourceROIStartX is odd
        {
          if (SourceROIWidth&1) //SourceROIWidth is odd
            {
              //Start point is odd; final point is even.
              do
                {
                  HrzCounter = DestROIWidth;
                  Value=(unsigned short)((unsigned int) pSource[0]+
                        (unsigned int) pSource[SourceWidthItems]);
                  *pDest = (unsigned char) (Value/2);
                  pSource++;
                  pDest++;
                  if (--HrzCounter)
                    {
                      do
                        {
                          Value=(unsigned short)((unsigned int) pSource[0]+(unsigned int) pSource[1]+
                                (unsigned int) pSource[SourceWidthItems]+(unsigned int) pSource[SourceWidthItems+1]);
                          *pDest = (unsigned char) (Value/4);
                          pSource += 2;
                          pDest++;
                        } while (--HrzCounter);
                    };
                  pSource+=AddToSource;
                  pDest  +=AddToDest;
                } while (--VrtCounter);
            }
          else //SourceROIWidth is even
            {
              //Start point is odd; final point is odd. Nonaligned loop.
              do
                {
                  HrzCounter = DestROIWidth;
                  do
                    {
                      Value=(unsigned short)((unsigned int) pSource[0]+(unsigned int) pSource[1]+
                            (unsigned int) pSource[SourceWidthItems]+(unsigned int) pSource[SourceWidthItems+1]);
                      *pDest = (unsigned char) (Value/4);
                      pSource += 2;
                      pDest++;
                    } while (--HrzCounter);
                  pSource+=AddToSource;
                  pDest  +=AddToDest;
                } while (--VrtCounter);
            };
        }
      else //SourceROIStartX is even
        {
          if (SourceROIWidth&1) //SourceROIWidth is odd
            {
              do
                {
                  HrzCounter = DestROIWidth;
                  if (--HrzCounter)
                    {
                      do
                        {
                          Value=(unsigned short)((unsigned int) pSource[0]+(unsigned int) pSource[1]+
                                (unsigned int) pSource[SourceWidthItems]+(unsigned int) pSource[SourceWidthItems+1]);
                          *pDest = (unsigned char) (Value/4);
                          pSource += 2;
                          pDest++;
                        } while (--HrzCounter);
                    };
                  Value=(unsigned short)((unsigned int) pSource[0]+
                        (unsigned int) pSource[SourceWidthItems]);
                  *pDest = (unsigned char) (Value/2);
                  pSource+=AddToSource+1;
                  pDest  +=AddToDest+1;
                } while (--VrtCounter);
            }
          else //SourceROIWidth is even
            {
              do
                {
                  HrzCounter = DestROIWidth;
                  do
                    {
                      Value=(unsigned short)((unsigned int) pSource[0]+(unsigned int) pSource[1]+
                            (unsigned int) pSource[SourceWidthItems]+(unsigned int) pSource[SourceWidthItems+1]);
                      *pDest = (unsigned char) (Value/4);
                      pSource += 2;
                      pDest++;
                    } while (--HrzCounter);
                  pSource+=AddToSource;
                  pDest  +=AddToDest;
                } while (--VrtCounter);
            };
        };
    };
}
