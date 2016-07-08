/*************************************************************
  FixFlt.c

  Project:     Computer Vision Library (CVL)
  
  Revision history:

  MM-YYYY    Author         Reason

  07-2000    SK,RE          initial C version
  10-2007    ANG            modified and updated C version
 *************************************************************/
#include <math.h>
#include "stddefs.h"

#ifdef _MSC_VER
  #pragma warning(disable:4146) // unary minus operator applied to unsigned type, result still unsigned
#endif

#define SET_VALUE(a,b)   \
{                       \
  tmpval = (b);         \
  if (tmpval<0)         \
    (a) = 0;            \
  else                  \
    if (tmpval>255)     \
      (a) = 255;        \
    else                \
      (a) = (unsigned char)tmpval; \
}
#define SET_VALUE_NEG(a,b)  \
  {                       \
  tmpval = (b);         \
  if (tmpval>0)         \
  (a) = 0;            \
    else                  \
    if (tmpval<(-255))  \
    (a) = 255;        \
      else                \
      (a) = (unsigned char)(-tmpval); \
  }

//#define CALC_NORM(v1,v2) ((abs(v1)+abs(v2))/2)
#define CALC_NORM(v1,v2) (sqrt((float)(v1*v1+v2*v2)))

#define SUM_TYPE  unsigned int
#define SET_VALUE_USHORT(a,b) (a) = (unsigned short)(b) 
#define CALC_NORM_USHORT(v1,v2) (sqrt((float)(v1*v1+v2*v2)))

// Function is separated into three parts called (UP), (MAIN), (LO).
// (UP) is realised if ROI touches the upper boundary of image, so the mask overflows
// image boundary. (LO) is for lower bound respectively. (MAIN) is for the inner part 
// of image, that is not influenced by boundaries.

void CVL_EdgeSobelFast_Byte(
  unsigned char*       Dest,
  const unsigned char* Source,
  const unsigned int   W,
  const unsigned int   H,
  const unsigned int   S,
  unsigned int         ROIStartX, 
  unsigned int         ROIStartY,
  unsigned int         ROIWidth,
  unsigned int         ROIHeight,
  int*                 MemoryBufferX, //SizeItems = 3*Width;
  int*                 MemoryBufferY) //SizeItems = 3*Width;
{
  int tmpval;
  // accumulator	       
  int* sum_x = MemoryBufferX;
  int* sum_y = MemoryBufferY;
  int *psum_y,*psum_x;
  int gx,gy;
  // pointers along y: source, destination, finishing source
  // pointers along x: source, destination, finishing source
  const unsigned char* src_point,*finish_x,*src_line,*finish_y;
  unsigned char *dst_line,*dst_point;
  size_t Width = W;
  size_t Height = H;
  size_t Stride = S;
  
  if (Dest&&Source&&sum_x&&sum_y&&(Width>2)&&(Height>2)&&ROIWidth&&ROIHeight&&
      ((ROIStartX+ROIWidth)<=Width)&&((ROIStartY+ROIHeight)<=Height))
  {
    // first, calculate borders where ROI touches image boundaries, if any
//-- TOP --------------------------------------------------------------------------------
    if (ROIStartY==0)
    { // ROI touches top of image
      src_point = Source;
      dst_point = Dest;
      if (ROIStartX+ROIWidth==Width)
      {
        // top-right corner
        gy =(int)src_point[Width       -2] +3*(int)src_point[Width       -1]
           -(int)src_point[Width+Stride-2] -3*(int)src_point[Width+Stride-1];
        gx = -3*(int)src_point[Width       -2] +3*(int)src_point[Width       -1]
               -(int)src_point[Width+Stride-2] +  (int)src_point[Width+Stride-1];

        SET_VALUE(dst_point[Width-1],(int)CALC_NORM(gy,gx));		
 		
        finish_x = src_point+Width-1;
      }
      else
        finish_x = src_point+ROIStartX+ROIWidth;

      if (ROIStartX==0)
      {
        // top-left corner of image
		    gy = 3*(int)src_point[  0  ] + 
                  +(int)src_point[  1  ] +
		        -3*(int)src_point[Stride  ]
			     - (int)src_point[Stride+1];
		    gx = -3*(int)src_point[  0  ] + 
                 +3*(int)src_point[  1  ] +
                  - (int)src_point[Stride  ]
                  + (int)src_point[Stride+1];

        SET_VALUE(dst_point[0],(int)CALC_NORM(gy,gx));
        src_point++;
        dst_point++;
      }
      else
      {
        src_point += ROIStartX;
        dst_point += ROIStartX;
      }
      // top line
      while (src_point<finish_x)
      {
        gy =    (int)src_point[      -1] +2*(int)src_point[       0] + (int)src_point[       1]
               -(int)src_point[Stride-1] -2*(int)src_point[Stride+0] - (int)src_point[Stride+1];
        gx = -3*(int)src_point[      -1] +3*(int)src_point[       1]
               -(int)src_point[Stride-1] +  (int)src_point[Stride+1];

        SET_VALUE(dst_point[0],(int)CALC_NORM(gy,gx));
        // move pointers
        src_point++;
        dst_point++;
      }
    }
//-- BOTTOM -----------------------------------------------------------------------------
    if (ROIStartY+ROIHeight==Height)
    { // ROI touches bottom of image
      src_point = Source+(Height-1)*Stride;
      dst_point = Dest+(Height-1)*Stride;
      if (ROIStartX+ROIWidth==Width)
      {
        // bottom-right corner
       gy =(int)src_point[Width-Stride-2] +3*(int)src_point[Width-Stride-1]
          -(int)src_point[Width       -2] -3*(int)src_point[Width       -1];
       gx = -(int)src_point[Width-Stride-2] +  (int)src_point[Width-Stride-1]
          -3*(int)src_point[Width       -2] +3*(int)src_point[Width       -1];

       SET_VALUE(dst_point[Width-1],(int)CALC_NORM(gy,gx));
       finish_x = src_point+Width-1;
      }
      else
        finish_x = src_point+ROIStartX+ROIWidth;
      if (ROIStartX==0)
      {
        // bottom-left corner of image
        gy = 3*(int)src_point[-Stride] +  (int)src_point[-Stride+1]
            -3*(int)src_point[      0] -  (int)src_point[        1];
        gx =  -(int)src_point[-Stride] +  (int)src_point[-Stride+1]
            -3*(int)src_point[      0] +3*(int)src_point[        1];
        SET_VALUE(dst_point[0],(int)CALC_NORM(gy,gx));
        src_point++;
        dst_point++;
      }
      else
      {
        src_point += ROIStartX;
        dst_point += ROIStartX;
      }
      // bottom line
      while (src_point<finish_x)
      {
        gy =  (int)src_point[-Stride-1] +2*(int)src_point[-Stride] + (int)src_point[-Stride+1]
             -(int)src_point[       -1] -2*(int)src_point[      0] - (int)src_point[        1];
        gx =  -(int)src_point[-Stride-1] +  (int)src_point[-Stride+1]
            -3*(int)src_point[       -1] +3*(int)src_point[        1];

        SET_VALUE(dst_point[0],(int)CALC_NORM(gy,gx));
        // move pointers
        src_point++;
        dst_point++;
      }
    }
//-- LEFT -------------------------------------------------------------------------------
    if (ROIStartX==0)
    { // ROI touches left border of image
      if (ROIStartY==0)
      { // exclude top-left corner of ROI - calculated already
        src_point = Source+Stride;
        dst_point = Dest+Stride;
      }
      else
      {
        src_point = Source+ROIStartY*Stride;
        dst_point = Dest+ROIStartY*Stride;
      }
      if (ROIStartY+ROIHeight==Height)
        // exclude bottom-left
        finish_y = Source+(Height-1)*Stride;
      else
        finish_y = Source+(ROIStartY+ROIHeight)*Stride;
      // left border line
      while (src_point<finish_y)
      {
        gy = 3*(int)src_point[-Stride ] +  (int)src_point[-Stride+1]
            -3*(int)src_point[ Stride ] -  (int)src_point[ Stride+1];
        gx =  -(int)src_point[-Stride ] +  (int)src_point[-Stride+1]
            -2*(int)src_point[       0] +2*(int)src_point[        1]
              -(int)src_point[ Stride ] +  (int)src_point[ Stride+1];

        SET_VALUE(dst_point[0],(int)CALC_NORM(gy,gx));
        // move pointers
        src_point += Stride;
        dst_point += Stride;
      }
    }
//-- RIGHT ------------------------------------------------------------------------------
    if (ROIStartX+ROIWidth==Width)
    { // ROI touches right border of image
      if (ROIStartY==0)
      { // exclude top-right corner of ROI - calculated already
        src_point = Source+Stride+Width-1;
        dst_point = Dest+Stride+Width-1;
      }
      else
      {
        src_point = Source+ROIStartY*Stride+Width-1;
        dst_point = Dest+ROIStartY*Stride+Width-1;
      }
      if (ROIStartY+ROIHeight==Height)
        // exclude bottom-left
        finish_y = Source+(Height-1)*Stride+Width-1;
      else
        finish_y = Source+(ROIStartY+ROIHeight)*Stride+Width-1;
      // right border line
      while (src_point<finish_y)
      {
          gy =    (int)src_point[-Stride-1] +3*(int)src_point[-Stride]
                 -(int)src_point[ Stride-1] -3*(int)src_point[ Stride];
          gx =   -(int)src_point[-Stride-1] +  (int)src_point[-Stride]
               -2*(int)src_point[       -1] +2*(int)src_point[      0]
                 -(int)src_point[ Stride-1] +  (int)src_point[ Stride];

          SET_VALUE(dst_point[0],(int)CALC_NORM(gy,gx));
        // move pointers
          src_point += Stride;
          dst_point += Stride;
      }
    }
//-- MAIN - not influenced by image borders ---------------------------------------------
    if (ROIStartY==0)   // touches top
    {
      ROIStartY = 1;
      ROIHeight--;
    }
    if (ROIStartY+ROIHeight==Height)  // touches bottom
      ROIHeight--;
    if (ROIStartX==0)   // touches left
    {
      ROIStartX = 1;
      ROIWidth--;
    }
    if (ROIStartX+ROIWidth==Width)  // touches right
      ROIWidth--;
    // cycle along y
    src_line = Source+ROIStartY*Stride+ROIStartX;
    dst_line = Dest+ROIStartY*Stride+ROIStartX;
    finish_y = src_line+ROIHeight*Stride;
    while (src_line!=finish_y)
    {
      src_point = src_line;
      dst_point = dst_line;
      finish_x  = src_line+ROIWidth;

      sum_y[0] =   (int)src_point [-Stride-1]                   
                  -(int)src_point [ Stride-1];
      sum_x[0] =  -(int)src_point [-Stride-1]
                -2*(int)src_point [       -1]
                  -(int)src_point [ Stride-1];	  
      sum_y[1] = 2*(int)src_point [-Stride]                 
                -2*(int)src_point [ Stride];
      sum_x[1]=   -(int)src_point [-Stride]
                -2*(int)src_point [      0]
                  -(int)src_point [ Stride];
      sum_y[2] =   (int)src_point [-Stride+1]+                            
                  -(int)src_point [ Stride+1];
      sum_x[2] =   (int)src_point [-Stride+1]+
                 2*(int)src_point [       +1]+
                   (int)src_point [ Stride+1];
      gy = sum_y[0] + sum_y[1] + sum_y[2];
      gx = sum_x[0] + sum_x[2];
    // write result
      SET_VALUE(dst_point[0],(int)CALC_NORM(gy,gx));        
      src_point++;
      dst_point++;
      psum_y = sum_y+3;
      psum_x = sum_x+3;
      while (src_point!=finish_x)
      {
        *psum_y     = *(psum_y-2)/2;
        *(psum_y+1) = *(psum_y-1)*2;
        *(psum_y+2) = (int)src_point[-Stride+1]+                            
                     -(int)src_point[ Stride+1] ;
        gy = *psum_y + *(psum_y+1) + *(psum_y+2);

        *psum_x     = *(psum_x-2);
        *(psum_x+1) = -*(psum_x-1);
        *(psum_x+2) = (int)src_point[-Stride+1]+
                    2*(int)src_point[       +1]+
                      (int)src_point[ Stride+1];
        gx = *psum_x + *(psum_x+2);

    //  write result
        SET_VALUE(dst_point[0],(int)CALC_NORM(gy,gx));
        psum_y += 3;
        psum_x += 3;        
        src_point++;
        dst_point++;
      }
      src_line += Stride;
      dst_line += Stride;
    }
  }
}

void CVL_EdgeSobelFast_ByteToUShort(
  unsigned short*    Dest,
  const unsigned char*  Source,
  const unsigned int W,
  const unsigned int H,
  const unsigned int S,
  unsigned int       ROIStartX, 
  unsigned int       ROIStartY,
  unsigned int       ROIWidth,
  unsigned int       ROIHeight,
  int*               MemoryBufferX, //SizeItems = 3*Width;
  int*               MemoryBufferY) //SizeItems = 3*Width;
{
  // accumulator	       
  int* sum_x = MemoryBufferX;
  int* sum_y = MemoryBufferY;
  int *psum_y,*psum_x;
  int gx,gy;
  size_t Width = W;
  size_t Height = H;
  size_t Stride = S;
  // pointers along y: source, destination, finishing source
  // pointers along x: source, destination, finishing source
  const unsigned char* src_point,*finish_x,*src_line,*finish_y;
  unsigned short *dst_line,*dst_point;
  if (Dest&&Source&&sum_x&&sum_y&&(Width>2)&&(Height>2)&&ROIWidth&&ROIHeight&&
      ((ROIStartX+ROIWidth)<=Width)&&((ROIStartY+ROIHeight)<=Height))
  {
    // first, calculate borders where ROI touches image boundaries, if any
//-- TOP --------------------------------------------------------------------------------
    if (ROIStartY==0)
    { // ROI touches top of image
      src_point = Source;
      dst_point = Dest;
      if (ROIStartX+ROIWidth==Width)
      {
        // top-right corner
        gy =(int)src_point[Width       -2] +3*(int)src_point[Width       -1]
           -(int)src_point[Width+Stride-2] -3*(int)src_point[Width+Stride-1];
        gx = -3*(int)src_point[Width       -2] +3*(int)src_point[Width       -1]
               -(int)src_point[Width+Stride-2] +  (int)src_point[Width+Stride-1];

        SET_VALUE_USHORT(dst_point[Width-1],(int)CALC_NORM_USHORT(gy,gx));		
 		
        finish_x = src_point+Width-1;
      }
      else
        finish_x = src_point+ROIStartX+ROIWidth;

      if (ROIStartX==0)
      {
        // top-left corner of image
		    gy = 3*(int)src_point[  0  ] + 
                  +(int)src_point[  1  ] +
		        -3*(int)src_point[Stride  ]
			     - (int)src_point[Stride+1];
		    gx = -3*(int)src_point[  0  ] + 
                 +3*(int)src_point[  1  ] +
                  - (int)src_point[Stride  ]
                  + (int)src_point[Stride+1];

        SET_VALUE_USHORT(dst_point[0],(int)CALC_NORM_USHORT(gy,gx));
        src_point++;
        dst_point++;
      }
      else
      {
        src_point += ROIStartX;
        dst_point += ROIStartX;
      }
      // top line
      while (src_point<finish_x)
      {
        gy =    (int)src_point[      -1] +2*(int)src_point[       0] + (int)src_point[       1]
               -(int)src_point[Stride-1] -2*(int)src_point[Stride+0] - (int)src_point[Stride+1];
        gx = -3*(int)src_point[      -1] +3*(int)src_point[       1]
               -(int)src_point[Stride-1] +  (int)src_point[Stride+1];

        SET_VALUE_USHORT(dst_point[0],(int)CALC_NORM_USHORT(gy,gx));
        // move pointers
        src_point++;
        dst_point++;
      }
    }
//-- BOTTOM -----------------------------------------------------------------------------
    if (ROIStartY+ROIHeight==Height)
    { // ROI touches bottom of image
      src_point = Source+(Height-1)*Stride;
      dst_point = Dest+(Height-1)*Stride;
      if (ROIStartX+ROIWidth==Width)
      {
        // bottom-right corner
       gy =(int)src_point[Width-Stride-2] +3*(int)src_point[Width-Stride-1]
          -(int)src_point[Width       -2] -3*(int)src_point[Width       -1];
//       gy = (int)(src_point[(size_t)Width-(size_t)Stride-2]);
//       gy += 3*(int)src_point[Width-Stride-1];
//       gy -= (int)src_point[Width       -2];
//       gy -= 3*(int)src_point[Width       -1];
       gx = -(int)src_point[Width-Stride-2] +  (int)src_point[Width-Stride-1]
          -3*(int)src_point[Width       -2] +3*(int)src_point[Width       -1];

       SET_VALUE_USHORT(dst_point[Width-1],(int)CALC_NORM_USHORT(gy,gx));
       finish_x = src_point+Width-1;
      }
      else
        finish_x = src_point+ROIStartX+ROIWidth;
      if (ROIStartX==0)
      {
        // bottom-left corner of image
        gy = 3*(int)src_point[-Stride] +  (int)src_point[-Stride+1]
            -3*(int)src_point[      0] -  (int)src_point[        1];
        gx =  -(int)src_point[-Stride] +  (int)src_point[-Stride+1]
            -3*(int)src_point[      0] +3*(int)src_point[        1];
        SET_VALUE_USHORT(dst_point[0],(int)CALC_NORM_USHORT(gy,gx));
        src_point++;
        dst_point++;
      }
      else
      {
        src_point += ROIStartX;
        dst_point += ROIStartX;
      }
      // bottom line
      while (src_point<finish_x)
      {
        gy =  (int)src_point[-Stride-1] +2*(int)src_point[-Stride] + (int)src_point[-Stride+1]
             -(int)src_point[       -1] -2*(int)src_point[      0] - (int)src_point[        1];
        gx =  -(int)src_point[-Stride-1] +  (int)src_point[-Stride+1]
            -3*(int)src_point[       -1] +3*(int)src_point[        1];

        SET_VALUE_USHORT(dst_point[0],(int)CALC_NORM_USHORT(gy,gx));
        // move pointers
        src_point++;
        dst_point++;
      }
    }
//-- LEFT -------------------------------------------------------------------------------
    if (ROIStartX==0)
    { // ROI touches left border of image
      if (ROIStartY==0)
      { // exclude top-left corner of ROI - calculated already
        src_point = Source+Stride;
        dst_point = Dest+Stride;
      }
      else
      {
        src_point = Source+ROIStartY*Stride;
        dst_point = Dest+ROIStartY*Stride;
      }
      if (ROIStartY+ROIHeight==Height)
        // exclude bottom-left
        finish_y = Source+(Height-1)*Stride;
      else
        finish_y = Source+(ROIStartY+ROIHeight)*Stride;
      // left border line
      while (src_point<finish_y)
      {
        gy = 3*(int)src_point[-Stride ] +  (int)src_point[-Stride+1]
            -3*(int)src_point[ Stride ] -  (int)src_point[ Stride+1];
        gx =  -(int)src_point[-Stride ] +  (int)src_point[-Stride+1]
            -2*(int)src_point[       0] +2*(int)src_point[        1]
              -(int)src_point[ Stride ] +  (int)src_point[ Stride+1];

        SET_VALUE_USHORT(dst_point[0],(int)CALC_NORM_USHORT(gy,gx));
        // move pointers
        src_point += Stride;
        dst_point += Stride;
      }
    }
//-- RIGHT ------------------------------------------------------------------------------
    if (ROIStartX+ROIWidth==Width)
    { // ROI touches right border of image
      if (ROIStartY==0)
      { // exclude top-right corner of ROI - calculated already
        src_point = Source+Stride+Width-1;
        dst_point = Dest+Stride+Width-1;
      }
      else
      {
        src_point = Source+ROIStartY*Stride+Width-1;
        dst_point = Dest+ROIStartY*Stride+Width-1;
      }
      if (ROIStartY+ROIHeight==Height)
        // exclude bottom-left
        finish_y = Source+(Height-1)*Stride+Width-1;
      else
        finish_y = Source+(ROIStartY+ROIHeight)*Stride+Width-1;
      // right border line
      while (src_point<finish_y)
      {
          gy =    (int)src_point[-Stride-1] +3*(int)src_point[-Stride]
                 -(int)src_point[ Stride-1] -3*(int)src_point[ Stride];
          gx =   -(int)src_point[-Stride-1] +  (int)src_point[-Stride]
               -2*(int)src_point[       -1] +2*(int)src_point[      0]
                 -(int)src_point[ Stride-1] +  (int)src_point[ Stride];

          SET_VALUE_USHORT(dst_point[0],(int)CALC_NORM_USHORT(gy,gx));
        // move pointers
          src_point += Stride;
          dst_point += Stride;
      }
    }
//-- MAIN - not influenced by image borders ---------------------------------------------
    if (ROIStartY==0)   // touches top
    {
      ROIStartY = 1;
      ROIHeight--;
    }
    if (ROIStartY+ROIHeight==Height)  // touches bottom
      ROIHeight--;
    if (ROIStartX==0)   // touches left
    {
      ROIStartX = 1;
      ROIWidth--;
    }
    if (ROIStartX+ROIWidth==Width)  // touches right
      ROIWidth--;
    // cycle along y
    src_line = Source+ROIStartY*Stride+ROIStartX;
    dst_line = Dest+ROIStartY*Stride+ROIStartX;
    finish_y = src_line+ROIHeight*Stride;
    while (src_line!=finish_y)
    {
      src_point = src_line;
      dst_point = dst_line;
      finish_x  = src_line+ROIWidth;

      sum_y[0] =   (int)src_point [-Stride-1]                   
                  -(int)src_point [ Stride-1];
      sum_x[0] =  -(int)src_point [-Stride-1]
                -2*(int)src_point [       -1]
                  -(int)src_point [ Stride-1];	  
      sum_y[1] = 2*(int)src_point [-Stride]                 
                -2*(int)src_point [ Stride];
      sum_x[1]=   -(int)src_point [-Stride]
                -2*(int)src_point [      0]
                  -(int)src_point [ Stride];
      sum_y[2] =   (int)src_point [-Stride+1]+                            
                  -(int)src_point [ Stride+1];
      sum_x[2] =   (int)src_point [-Stride+1]+
                 2*(int)src_point [       +1]+
                   (int)src_point [ Stride+1];
      gy = sum_y[0] + sum_y[1] + sum_y[2];
      gx = sum_x[0] + sum_x[2];
    // write result
      SET_VALUE_USHORT(dst_point[0],(int)CALC_NORM_USHORT(gy,gx));        
      src_point++;
      dst_point++;
      psum_y = sum_y+3;
      psum_x = sum_x+3;
      while (src_point!=finish_x)
      {
        *psum_y     = *(psum_y-2)/2;
        *(psum_y+1) = *(psum_y-1)*2;
        *(psum_y+2) = (int)src_point[-Stride+1]+                            
                     -(int)src_point[ Stride+1] ;
        gy = *psum_y + *(psum_y+1) + *(psum_y+2);

        *psum_x     = *(psum_x-2);
        *(psum_x+1) = -*(psum_x-1);
        *(psum_x+2) = (int)src_point[-Stride+1]+
                    2*(int)src_point[       +1]+
                      (int)src_point[ Stride+1];
        gx = *psum_x + *(psum_x+2);

    //  write result
        SET_VALUE_USHORT(dst_point[0],(int)CALC_NORM_USHORT(gy,gx));
        psum_y += 3;
        psum_x += 3;        
        src_point++;
        dst_point++;
      }
      src_line += Stride;
      dst_line += Stride;
    }
  }
}
