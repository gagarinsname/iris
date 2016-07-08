/**********************************************************************    
  DrawCirc.c
 
   Project
      Computer Vision Library (CVL)

   Author
      ANG,SVY

 **********************************************************************/

// Brezenham's algorithm
void CVL_FillCircleInRectUShort(
  unsigned short* Source,
  unsigned int   Width,
  unsigned int   Height,
  unsigned int   WidthItems,
  unsigned int   ROIStartX,
  unsigned int   ROIStartY,
  unsigned int   ROIWidth,
  unsigned int   ROIHeight,
  int            XCenter,
  int            YCenter,
  unsigned int   Radius,
  unsigned short ValueToSet)
{
  unsigned short* pSource;
  unsigned short* pTopSource;
  unsigned short* pBottomSource;
  unsigned short* pTopSource2;
  unsigned short* pBottomSource2;
//  unsigned int   TotalPointsValue=0;
//  unsigned int   PoinsMarkedValue=0;
  unsigned int   X;
  unsigned int   Y;
  int            LeftROIBound =((int ) ROIStartX)-XCenter;
  int            RightROIBound=((int ) (ROIStartX+ROIWidth))-XCenter;
  int            TopROIBound =((int ) ROIStartY)-YCenter;
  int            BottomROIBound=((int ) (ROIStartY+ROIHeight))-YCenter;
  int            P;
  int            XLeft;
  int            XRight;
  unsigned int   BytesToSet;

  if (!Width || !Height ||(WidthItems<Width))
    return;

  if (!Radius)
    Radius++;    
   X = Radius;
   Y = 0;
   P = 3 - (int ) (Radius<<1);
   pTopSource    =Source+WidthItems*YCenter+XCenter;
   pBottomSource =pTopSource;
   pTopSource2   =Source+WidthItems*(YCenter-(int ) Radius)+XCenter;
   pBottomSource2=Source+WidthItems*(YCenter+(int ) Radius)+XCenter;
   if (Y<X)
     {
       do
       {
         if (((-(int ) Y)>=TopROIBound)&&((-(int ) Y)<(int ) BottomROIBound))
           {
             XLeft=(-(int ) X);
             if (XLeft<LeftROIBound)
               XLeft=LeftROIBound;
             XRight=((int ) X);
             if (XRight>RightROIBound)
               XRight=RightROIBound;
             if (XLeft<=XRight)
               {
                 BytesToSet=(unsigned int ) (XRight-XLeft)+1;
                 pSource=pTopSource+XLeft-1;
                 do
                   {
                     pSource[BytesToSet]=ValueToSet;
                   } while (--BytesToSet);
               };
           };
         if (((int ) Y>=TopROIBound)&&((int ) Y<(int ) BottomROIBound))
           {
             XLeft=(-(int ) X);
             if (XLeft<LeftROIBound)
               XLeft=LeftROIBound;
             XRight=((int ) X);
             if (XRight>RightROIBound)
               XRight=RightROIBound;
             if (XLeft<=XRight)
               {
                 BytesToSet=(unsigned int ) (XRight-XLeft)+1;
                 pSource=pBottomSource+XLeft-1;
                 do
                   {
                     pSource[BytesToSet]=ValueToSet;
                   } while (--BytesToSet);
               };
           };
         if (((-(int ) X)>=TopROIBound)&&((-(int ) X)<(int ) BottomROIBound))
           {
             XLeft=(-(int ) Y);
             if (XLeft<LeftROIBound)
               XLeft=LeftROIBound;
             XRight=((int ) Y);
             if (XRight>RightROIBound)
               XRight=RightROIBound;
             if (XLeft<=XRight)
               {
                 BytesToSet=(unsigned int ) (XRight-XLeft)+1;
                 pSource=pTopSource2+XLeft-1;
                 do
                   {
                     pSource[BytesToSet]=ValueToSet;
                   } while (--BytesToSet);
               };
           };
         if (((int ) X>=TopROIBound)&&((int ) X<(int ) BottomROIBound))
           {
             XLeft=(-(int ) Y);
             if (XLeft<LeftROIBound)
               XLeft=LeftROIBound;
             XRight=((int ) Y);
             if (XRight>RightROIBound)
               XRight=RightROIBound;
             if (XLeft<=XRight)
               {
                 BytesToSet=(unsigned int ) (XRight-XLeft)+1;
                 pSource=pBottomSource2+XLeft-1;
                 do
                   {
                     pSource[BytesToSet]=ValueToSet;
                   } while (--BytesToSet);
               };
           };
         if (P<0)
           {
             P+=4*((int ) Y)+6;
           }
         else
           {
             P+=4*((int ) Y-(int ) X)+10;
             X--;
             pTopSource2   +=WidthItems;
             pBottomSource2-=WidthItems;
           };
         pTopSource   -=WidthItems;
         pBottomSource+=WidthItems;
       } while (++Y<X);
   };
   if (X==Y)
     {
       if (((-(int ) Y)>=TopROIBound)&&((-(int ) Y)<(int ) BottomROIBound))
         {
           XLeft=(-(int ) X);
           if (XLeft<LeftROIBound)
             XLeft=LeftROIBound;
           XRight=((int ) X);
           if (XRight>RightROIBound)
             XRight=RightROIBound;
           if (XLeft<=XRight)
             {
               BytesToSet=(unsigned int ) (XRight-XLeft)+1;
               pSource=pTopSource+XLeft-1;
               do
                 {
                   pSource[BytesToSet]=ValueToSet;
                 } while (--BytesToSet);
             };
         };
       if (((int ) Y>=TopROIBound)&&((int ) Y<(int ) BottomROIBound))
         {
           XLeft=(-(int ) X);
           if (XLeft<LeftROIBound)
             XLeft=LeftROIBound;
           XRight=((int ) X);
           if (XRight>RightROIBound)
             XRight=RightROIBound;
           if (XLeft<=XRight)
             {
               BytesToSet=(unsigned int ) (XRight-XLeft)+1;
               pSource=pBottomSource+XLeft-1;
               do
                 {
                   pSource[BytesToSet]=ValueToSet;
                 } while (--BytesToSet);
             };
         };
     };
   
}

//Brezenhame's algorithm
void CVL_FillHorHalfCircleRingInRectUShort(
  unsigned short* Source,
  unsigned int   Width,
  unsigned int   Height,
  unsigned int   WidthItems,
  unsigned int   ROIStartX,
  unsigned int   ROIStartY,
  unsigned int   ROIWidth,
  unsigned int   ROIHeight,
  int            XCenter,
  int            YCenter,
  unsigned int   StartRadius,
  unsigned int   FinalRadius,
  unsigned short ValueToSet)
{
  unsigned short* pSource;
  unsigned short* pTopSource;
  unsigned short* pBottomSource;
  unsigned int   X,Xs;
  unsigned int   Y;
  int            LeftROIBound =((int ) ROIStartX)-XCenter;
  int            RightROIBound=((int ) (ROIStartX+ROIWidth))-XCenter;
  int            TopROIBound =((int ) ROIStartY)-YCenter;
  int            BottomROIBound=((int ) (ROIStartY+ROIHeight))-YCenter;
  int            P,Ps;
  int            XLeft;
  int            XRight;
  unsigned int   BytesToSet;

  if (!Width || !Height ||(WidthItems<Width))
    return;

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
      
      //if (Y<X)
        {
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
                      BytesToSet=(unsigned int ) (XRight-XLeft)+1;
                      pSource=pTopSource+XLeft-1;
                      do
                        {
                          pSource[BytesToSet]=ValueToSet;
                        } while (--BytesToSet);
                  };
                  XLeft=((int ) Xs);
                  if (XLeft<LeftROIBound)
                    XLeft=LeftROIBound;
                  XRight=((int ) X);
                  if (XRight>RightROIBound)
                    XRight=RightROIBound;
                  if ((XLeft<=RightROIBound)&&(XLeft<=XRight))
                    {
                      BytesToSet=(unsigned int ) (XRight-XLeft)+1;
                      pSource=pTopSource+XLeft-1;
                      do
                        {
                          pSource[BytesToSet]=ValueToSet;
                        } while (--BytesToSet);
                    };
                };
                             
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
                      BytesToSet=(unsigned int ) (XRight-XLeft)+1;
                      pSource=pBottomSource+XLeft-1;
                      do
                        {
                          pSource[BytesToSet]=ValueToSet;
                        } while (--BytesToSet);
                  };
                  XLeft=((int ) Xs);
                  if (XLeft<LeftROIBound)
                    XLeft=LeftROIBound;
                  XRight=((int ) X);
                  if (XRight>RightROIBound)
                    XRight=RightROIBound;
                  if ((XLeft<=RightROIBound)&&(XLeft<=XRight))
                    {
                      BytesToSet=(unsigned int ) (XRight-XLeft)+1;
                      pSource=pBottomSource+XLeft-1;
                      do
                        {
                          pSource[BytesToSet]=ValueToSet;
                        } while (--BytesToSet);
                    };                     
                };                                 

              if (P<0)
              {
                 P+=4*((int ) Y)+6;
              }
              else
              {
                 P+=4*((int ) Y-(int ) X)+10;
                 X--;                     
              };
              
              if (StartRadius)
              {
               if (Ps<0)
               {
                  Ps+=4*((int ) Y)+6;                     
               }
               else
               {
                  Ps+=4*((int ) Y-(int ) Xs)+10;
                  Xs--;                      
               };                  
              };

              pTopSource   -=WidthItems;
              pBottomSource+=WidthItems;               
            } while (++Y<=X);
        };
      /*
      if (X==Y)
        {
          if (((-(int ) Y)>=TopROIBound)&&((-(int ) Y)<(int ) BottomROIBound))
            {
              XLeft=(-(int ) X);
              if (XLeft<LeftROIBound)
                XLeft=LeftROIBound;
              XRight=((int ) X);
              if (XRight>RightROIBound)
                XRight=RightROIBound;
              if (XLeft<=XRight)
                {
                  BytesToSet=(unsigned int ) (XRight-XLeft)+1;
                  pSource=pTopSource+XLeft-1;
                  do
                    {
                      pSource[BytesToSet]=ValueToSet;
                    } while (--BytesToSet);
                };
            };
          if (((int ) Y>=TopROIBound)&&((int ) Y<(int ) BottomROIBound))
            {
              XLeft=(-(int ) X);
              if (XLeft<LeftROIBound)
                XLeft=LeftROIBound;
              XRight=((int ) X);
              if (XRight>RightROIBound)
                XRight=RightROIBound;
              if (XLeft<=XRight)
                {
                  BytesToSet=(unsigned int ) (XRight-XLeft)+1;
                  pSource=pBottomSource+XLeft-1;
                  do
                    {
                      pSource[BytesToSet]=ValueToSet;
                    } while (--BytesToSet);
                };
            };
        };
        */
     
    };
}
