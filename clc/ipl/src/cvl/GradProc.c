/**********************************************************************    
  GradProc.c
 
   Project
      Computer Vision Library (CVL)

   Author
      Alexander N. Gneushev

  Description
      Gradient functions

 **********************************************************************/

#include <math.h>
#include "stddefs.h"

#ifdef _MSC_VER
  #pragma warning(disable:4146) // unary minus operator applied to unsigned type, result still unsigned
#endif

#define LOOP_SUM_GARD_IN_HOR_RING     \
do                                    \
{                                     \
 pSrc = pSource+BytesToSet;           \
 SourceValue = pSrc[0];               \
 if (SourceValue)                     \
{                                     \
  GradArea += SourceValue;            \
  if (((SourceValue > pSrc[1])&&(SourceValue > pSrc[-1]))||                       \
      ((SourceValue > pSrc[-WidthItems+1])&&(SourceValue > pSrc[WidthItems-1]))|| \
      ((SourceValue > pSrc[-WidthItems])&&(SourceValue > pSrc[WidthItems]))||     \
      ((SourceValue > pSrc[WidthItems+1])&&(SourceValue > pSrc[-WidthItems-1])))  \
  {                                                                               \
   MaxGradArea += SourceValue;                                                    \
  }                                                                               \
 }                                                                                \
} while (--BytesToSet)

void CVL_CalcMaxGradAreaUShort(
  unsigned int*         pMaxGradArea,
  unsigned int*         pGradArea,
  const unsigned short* Source,                   
  unsigned int          ROIStartX,
  unsigned int          ROIStartY,
  unsigned int          ROIWidth,
  unsigned int          ROIHeight,
  unsigned int          WI)
{
 unsigned int         i,j;
 const unsigned short* pSource  = Source  + ROIStartX+WI*ROIStartY;
 const unsigned int   AddToPtr = WI-ROIWidth;
 unsigned short       SourceValue;
 unsigned int         MaxGradArea = 0;
 unsigned int         GradArea    = 0;
 size_t WidthItems = WI;


 if (pMaxGradArea)
   *pMaxGradArea = 0;

 if (pGradArea)
   *pGradArea = 0;

 if (pMaxGradArea&&Source&&ROIWidth&&ROIHeight&&
     (ROIStartX+ROIWidth<=WidthItems))
 {
//////////////////////////    
     SourceValue = pSource[0];
     if (SourceValue)
     {
       GradArea += SourceValue;
        if (((SourceValue > pSource[1]))||
            ((SourceValue > pSource[WidthItems]))||
            ((SourceValue > pSource[WidthItems+1])))
                 MaxGradArea += SourceValue;
     }    
     pSource++;
    ///
    i = ROIWidth-2;
    do
    {    
      SourceValue = pSource[0];      
      if (SourceValue)
      {
        GradArea += SourceValue;
        if (((SourceValue > pSource[1])&&(SourceValue > pSource[-1]))||
            ((SourceValue > pSource[WidthItems-1]))||
            ((SourceValue > pSource[WidthItems]))||
            ((SourceValue > pSource[WidthItems+1])))
                 MaxGradArea += SourceValue;
      }    
      pSource++;
    }while (--i);
    ///  
     SourceValue = pSource[0];     
     if (SourceValue)
     {
       GradArea += SourceValue;
        if (((SourceValue > pSource[-1]))||
            ((SourceValue > pSource[WidthItems-1]))||
            ((SourceValue > pSource[WidthItems])))
                       MaxGradArea += SourceValue;
     }   
     pSource++;
    ///   
    pSource  += AddToPtr;

//////////////////////////
   j = ROIHeight-2;
   do
   {
     ///     
     SourceValue = pSource[0];
     if (SourceValue)
     {
       GradArea += SourceValue;
        if (((SourceValue > pSource[1]))||
            ((SourceValue > pSource[-WidthItems+1]))||
            ((SourceValue > pSource[-WidthItems])&&(SourceValue > pSource[WidthItems]))||
            ((SourceValue > pSource[WidthItems+1])))
                 MaxGradArea += SourceValue;
     }    
     pSource++;
    ///
    i = ROIWidth-2;
    do
    {      
      SourceValue = pSource[0];
      if (SourceValue)
      {
        GradArea += SourceValue;
        if (((SourceValue > pSource[1])&&(SourceValue > pSource[-1]))||
            ((SourceValue > pSource[-WidthItems+1])&&(SourceValue > pSource[WidthItems-1]))||
            ((SourceValue > pSource[-WidthItems])&&(SourceValue > pSource[WidthItems]))||
            ((SourceValue > pSource[WidthItems+1])&&(SourceValue > pSource[-WidthItems-1])))
                 MaxGradArea += SourceValue;
      }     
      pSource++;
    }while (--i);
    ///    
     SourceValue = pSource[0];
     if (SourceValue)
     {
       GradArea += SourceValue;
        if (((SourceValue > pSource[-1]))||
            ((SourceValue > pSource[WidthItems-1]))||
            ((SourceValue > pSource[-WidthItems])&&(SourceValue > pSource[WidthItems]))||
            ((SourceValue > pSource[-WidthItems-1])))
                       MaxGradArea += SourceValue;
     }    
     pSource++;
    ///
    pSource  += AddToPtr;
   }while(--j);
//////////////////////////     
     SourceValue = pSource[0];
     if (SourceValue)
     {
       GradArea += SourceValue;
        if (((SourceValue > pSource[1]))||
            ((SourceValue > pSource[-WidthItems+1]))||
            ((SourceValue > pSource[-WidthItems])))
                 MaxGradArea += SourceValue;
     }    
     pSource++;
    ///
    i = ROIWidth-2;
    do
    {      
      SourceValue = pSource[0];
      if (SourceValue)
      {
        GradArea += SourceValue;
        if (((SourceValue > pSource[1])&&(SourceValue > pSource[-1]))||
            ((SourceValue > pSource[-WidthItems+1]))||
            ((SourceValue > pSource[-WidthItems]))||
            ((SourceValue > pSource[-WidthItems-1])))
                 MaxGradArea += SourceValue;
      }     
      pSource++;
    }while (--i);
    ///     
     SourceValue = pSource[0];
     if (SourceValue)
     {
       GradArea += SourceValue;
        if (((SourceValue > pSource[-1]))||
            ((SourceValue > pSource[-WidthItems]))||
            ((SourceValue > pSource[-WidthItems-1])))
                       MaxGradArea += SourceValue;
     }     

     *pMaxGradArea = MaxGradArea;
     if (pGradArea)
       *pGradArea = GradArea;
 }
}

void CVL_CalcMaxGradAreaInCircHorHalfRingUShort(
  unsigned int*         pMaxGradArea,
  unsigned int*         pGradArea,
  const unsigned short* Source,                   
  unsigned int          ROIStartX,
  unsigned int          ROIStartY,
  unsigned int          ROIWidth,
  unsigned int          ROIHeight,
  unsigned int          WI,
  int                   XCenter,
  int                   YCenter,
  unsigned int          StartRadius,
  unsigned int          FinalRadius)
{
  const unsigned short* pSource;
  const unsigned short* pTopSource;
  const unsigned short* pBottomSource;
  const unsigned short* pSrc;
  unsigned int   X,Xs;
  unsigned int   Y;
  int            LeftROIBound =((int ) ROIStartX+1)-XCenter;
  int            RightROIBound=((int ) (ROIStartX+ROIWidth-1))-XCenter;
  int            TopROIBound =((int ) ROIStartY+1)-YCenter;
  int            BottomROIBound=((int ) (ROIStartY+ROIHeight-1))-YCenter;
  int            P,Ps;
  int            XLeft;
  int            XRight;
  unsigned int   BytesToSet;
  unsigned int   GradArea = 0;
  unsigned int   MaxGradArea = 0;
  unsigned short SourceValue;
  size_t          WidthItems = WI;

    if (!FinalRadius)
      FinalRadius++;
    if (StartRadius<=FinalRadius)
    {
      Y = 0;
      X = FinalRadius;
      P = 3 - (int ) (FinalRadius<<1);
      Xs = StartRadius;
      Ps = 3 - (int ) (StartRadius<<1);
      pTopSource    =Source+WidthItems*YCenter+XCenter;
      pBottomSource =pTopSource;      
      do
      {
        if (StartRadius&&(Y >= Xs))
          Xs = Y;
        if (((-(int ) Y)>=TopROIBound)&&((-(int ) Y)<(int ) BottomROIBound))
        {
          XLeft=(-(int ) X);
          if (XLeft<LeftROIBound)
            XLeft=LeftROIBound;
          XRight = (-(int ) Xs);
          if (XRight>RightROIBound)
            XRight=RightROIBound;
          if ((XRight>=LeftROIBound)&&(XLeft<=XRight))
          {
            BytesToSet=(unsigned int) (XRight-XLeft)+1;
            pSource=pTopSource+XLeft-1;
            LOOP_SUM_GARD_IN_HOR_RING;
          }
          XLeft=((int ) Xs);
          if (XLeft<LeftROIBound)
            XLeft=LeftROIBound;
          XRight=((int ) X);
          if (XRight>RightROIBound)
            XRight=RightROIBound;
          if ((XLeft<=RightROIBound)&&(XLeft<=XRight))
          {
            BytesToSet=(unsigned int) (XRight-XLeft)+1;
            pSource=pTopSource+XLeft-1;
            LOOP_SUM_GARD_IN_HOR_RING;
          }
        }
        if (((int ) Y>=TopROIBound)&&((int ) Y<(int ) BottomROIBound))
        {
          XLeft=(-(int ) X);
          if (XLeft<LeftROIBound)
            XLeft=LeftROIBound;
          XRight = (-(int ) Xs);
          if (XRight>RightROIBound)
            XRight=RightROIBound;
          if ((XRight>=LeftROIBound)&&(XLeft<=XRight))
          {
            BytesToSet=(unsigned int) (XRight-XLeft)+1;
            pSource=pBottomSource+XLeft-1;
            LOOP_SUM_GARD_IN_HOR_RING;
          }
          XLeft=((int ) Xs);
          if (XLeft<LeftROIBound)
            XLeft=LeftROIBound;
          XRight=((int ) X);
          if (XRight>RightROIBound)
            XRight=RightROIBound;
          if ((XLeft<=RightROIBound)&&(XLeft<=XRight))
          {
            BytesToSet=(unsigned int) (XRight-XLeft)+1;
            pSource=pBottomSource+XLeft-1;
            LOOP_SUM_GARD_IN_HOR_RING;
          }
        }
        if (P<0)
          P+=4*((int ) Y)+6;
        else
        {
          P+=4*((int ) Y-(int ) X)+10;
          X--;                     
        }
        if (StartRadius)
        {
          if (Ps<0)
            Ps+=4*((int ) Y)+6;                     
          else
          {
            Ps+=4*((int ) Y-(int ) Xs)+10;
            Xs--;                      
          }
        }
        pTopSource   -=WidthItems;
        pBottomSource+=WidthItems;               
      }
      while (++Y<=X);
    }
    if (pMaxGradArea)
      *pMaxGradArea = MaxGradArea;
    if (pGradArea)
      *pGradArea = GradArea;
}
