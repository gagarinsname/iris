/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:     13 April 2001                                       */
/*      Revision:    1.0.00                                              */
/*      Purpose:     rotate and resize a ROI of image                    */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#define __FILENUM__ 4 // __FILENUM__TAG4

#include <math.h>
#include "bpl.h"
#include "ipl.h"

// rotate and resize ROI with bi-linear interpolation
MIA_RESULT_CODE IPL_GEOM_RotateResizeBLI(
  unsigned char* dst,       // destination image
  int dststr,               // stride in bytes
  const unsigned char* src, // source image
  int srcstr,               // stride in bytes
  int xs,                   // source image size
  int ys,                   //
  float x,                  // source ROI
  float y,                  //
  float w,                  //
  float h,                  //
  float a,                  // angle of source ROI tilt
  int W,                    // destination ROI size
  int H)                    //
{
  /*register*/ int grossx,grossy;
  int xdX,ydX,xdY,ydY,i,j,add2dst,sum;
  const unsigned char* ptr;

    // +++ here one should insert the separating and different treatment
    //     of three cases:
    //      1) no rotation, no resize - simple copying ROI
    //      2) no rotation - simple resize
    //      3) rotation is present - this is the most universal most complex
    //         and thus slowest method. It is implemented now.
    //         Furthermore, rotations by degrees of pi/2 without resising
    //         can be outgrouped to one more case - mirroring.
    // ---

    add2dst = dststr-W;
    // form deltas
    xdX =  (int)((MIA_FIX20_cos((int)(a*(1<<20)))*w)/(W<<4));
    ydX =  (int)((MIA_FIX20_sin((int)(a*(1<<20)))*w)/(W<<4));
    xdY = -(int)((MIA_FIX20_sin((int)(a*(1<<20)))*h)/(H<<4));
    ydY =  (int)((MIA_FIX20_cos((int)(a*(1<<20)))*h)/(H<<4));
    // main cycle
    for (i=0;i<H;i++)
    {
      grossx = xdY*i+((int)(x*0x00010000));
      grossy = ydY*i+((int)(y*0x00010000));
      for (j=0;j<W;j++)
      {
        if ((grossx>=0)&&
            ((grossx>>16)+1<(int)xs)&&
            (grossy>=0)&&
            ((grossy>>16)+1<(int)ys) )
        {
          sum = 0x8000;
          ptr = src+ ((grossy>>16)*srcstr)+(grossx>>16);
          sum += ptr[0]*        ((0x10000-(grossx&0xffff))>>8)*
                                ((0x10000-(grossy&0xffff))>>8);
          sum += ptr[1]*        ((        (grossx&0xffff))>>8)*
                                ((0x10000-(grossy&0xffff))>>8);
          sum += ptr[srcstr]*   ((0x10000-(grossx&0xffff))>>8)*
                                ((        (grossy&0xffff))>>8);
          sum += ptr[srcstr+1]* ((        (grossx&0xffff))>>8)*
                                ((        (grossy&0xffff))>>8);
          *dst = (unsigned char)(sum>>16);
        }
        else
          *dst = 0;
        dst++;
        grossx += xdX;
        grossy += ydX;
      }
      dst += add2dst;
    }
    return ERR_OK;
}

// rotate and resize ROI with bi-linear interpolation
MIA_RESULT_CODE IPL_GEOM_AffineTransformBLI_RGB(
  unsigned char* dst,         // destination image
  unsigned int dststr,        // stride in bytes
  const unsigned char* src,   // source image
  unsigned int srcstr,        // stride in bytes
  unsigned int xs,            // source image size
  unsigned int ys,            //
  double x,                   // source ROI
  double y,                   //
  const double* A,            // affine transformation matrix
  unsigned int W,             // destination ROI size
  unsigned int H)             //
{
  /*register*/ int grossx,grossy;
  int xdX,ydX,xdY,ydY;
  unsigned int sumX,i,j,add2dst;
  const unsigned char* ptr;

    // check arguments
    if ((src==NULL)||(dst==NULL)||(A==NULL))
      return ERR_GEN_NULLPOINTER;
    add2dst = dststr-W*3;
    // form deltas
    xdX = (int)(A[0]*0x00010000+.5);// (int)((cos(a)*w*0x00010000)/W+.5);
    ydX = (int)(A[2]*0x00010000+.5);// (int)((sin(a)*w*0x00010000)/W+.5);
    xdY = (int)(A[1]*0x00010000+.5);//-(int)((sin(a)*h*0x00010000)/H+.5);
    ydY = (int)(A[3]*0x00010000+.5);// (int)((cos(a)*h*0x00010000)/H+.5);
    // main cycle
    for (i=0;i<H;i++)
    {
      grossx = xdY*i+((int)(x*0x00010000));
      grossy = ydY*i+((int)(y*0x00010000));
      for (j=0;j<W;j++)
      {
        if ((grossx>=0)&&
            ((grossx>>16)+1<(int)xs)&&
            (grossy>=0)&&
            ((grossy>>16)+1<(int)ys) )
        {
          ptr = src+ ((grossy>>16)*srcstr)+(grossx>>16)*3;
          // R
          sumX = 0x8000;
          sumX += ptr[0]*       ((0x10000-(grossx&0xffff))>>8)*
                                ((0x10000-(grossy&0xffff))>>8);
          sumX += ptr[3]*       ((        (grossx&0xffff))>>8)*
                                ((0x10000-(grossy&0xffff))>>8);
          sumX += ptr[srcstr+0]*((0x10000-(grossx&0xffff))>>8)*
                                ((        (grossy&0xffff))>>8);
          sumX += ptr[srcstr+3]*((        (grossx&0xffff))>>8)*
                                ((        (grossy&0xffff))>>8);
          dst[0] = (unsigned char)(sumX>>16);
          // G
          sumX = 0x8000;
          sumX += ptr[1]*       ((0x10000-(grossx&0xffff))>>8)*
                                ((0x10000-(grossy&0xffff))>>8);
          sumX += ptr[4]*       ((        (grossx&0xffff))>>8)*
                                ((0x10000-(grossy&0xffff))>>8);
          sumX += ptr[srcstr+1]*((0x10000-(grossx&0xffff))>>8)*
                                ((        (grossy&0xffff))>>8);
          sumX += ptr[srcstr+4]*((        (grossx&0xffff))>>8)*
                                ((        (grossy&0xffff))>>8);
          dst[1] = (unsigned char)(sumX>>16);
          // B
          sumX = 0x8000;
          sumX += ptr[2]*       ((0x10000-(grossx&0xffff))>>8)*
                                ((0x10000-(grossy&0xffff))>>8);
          sumX += ptr[5]*       ((        (grossx&0xffff))>>8)*
                                ((0x10000-(grossy&0xffff))>>8);
          sumX += ptr[srcstr+2]*((0x10000-(grossx&0xffff))>>8)*
                                ((        (grossy&0xffff))>>8);
          sumX += ptr[srcstr+5]*((        (grossx&0xffff))>>8)*
                                ((        (grossy&0xffff))>>8);
          dst[2] = (unsigned char)(sumX>>16);
        }
        else
          dst[0] = dst[1] = dst[2] = 0;
        dst += 3;
        grossx += xdX;
        grossy += ydX;
      }
      dst += add2dst;
    }
    return ERR_OK;
}

// rotate and resize ROI with bi-linear interpolation
MIA_RESULT_CODE IPL_GEOM_AffineTransformBLI_RGBtoGRAY(
  unsigned char* dst,         // destination image
  unsigned int dststr,        // stride in bytes
  const unsigned char* src,   // source image
  unsigned int srcstr,        // stride in bytes
  unsigned int xs,            // source image size
  unsigned int ys,            //
  double x,                   // source ROI
  double y,                   //
  const double* A,            // affine transformation matrix
  unsigned int W,             // destination ROI size
  unsigned int H)             //
{
  /*register*/ int grossx,grossy;
  int xdX,ydX,xdY,ydY;
  unsigned int sumX,i,j,add2dst;
  const unsigned char* ptr;

    // check arguments
    if ((src==NULL)||(dst==NULL)||(A==NULL))
      return ERR_GEN_NULLPOINTER;
    add2dst = dststr-W;
    // form deltas
    xdX = (int)(A[0]*0x00010000+.5);// (int)((cos(a)*w*0x00010000)/W+.5);
    ydX = (int)(A[2]*0x00010000+.5);// (int)((sin(a)*w*0x00010000)/W+.5);
    xdY = (int)(A[1]*0x00010000+.5);//-(int)((sin(a)*h*0x00010000)/H+.5);
    ydY = (int)(A[3]*0x00010000+.5);// (int)((cos(a)*h*0x00010000)/H+.5);
    // main cycle
    for (i=0;i<H;i++)
    {
      grossx = xdY*i+((int)(x*0x00010000));
      grossy = ydY*i+((int)(y*0x00010000));
      for (j=0;j<W;j++)
      {
        if ((grossx>=0)&&
            ((grossx>>16)+1<(int)xs)&&
            (grossy>=0)&&
            ((grossy>>16)+1<(int)ys) )
        {
          ptr = src+ ((grossy>>16)*srcstr)+(grossx>>16)*3;
          sumX = 0x8000*3;
          // R
          sumX += ptr[0]*       ((0x10000-(grossx&0xffff))>>8)*
                                ((0x10000-(grossy&0xffff))>>8);
          sumX += ptr[3]*       ((        (grossx&0xffff))>>8)*
                                ((0x10000-(grossy&0xffff))>>8);
          sumX += ptr[srcstr+0]*((0x10000-(grossx&0xffff))>>8)*
                                ((        (grossy&0xffff))>>8);
          sumX += ptr[srcstr+3]*((        (grossx&0xffff))>>8)*
                                ((        (grossy&0xffff))>>8);
          // G
          sumX += ptr[1]*       ((0x10000-(grossx&0xffff))>>8)*
                                ((0x10000-(grossy&0xffff))>>8);
          sumX += ptr[4]*       ((        (grossx&0xffff))>>8)*
                                ((0x10000-(grossy&0xffff))>>8);
          sumX += ptr[srcstr+1]*((0x10000-(grossx&0xffff))>>8)*
                                ((        (grossy&0xffff))>>8);
          sumX += ptr[srcstr+4]*((        (grossx&0xffff))>>8)*
                                ((        (grossy&0xffff))>>8);
          // B
          sumX += ptr[2]*       ((0x10000-(grossx&0xffff))>>8)*
                                ((0x10000-(grossy&0xffff))>>8);
          sumX += ptr[5]*       ((        (grossx&0xffff))>>8)*
                                ((0x10000-(grossy&0xffff))>>8);
          sumX += ptr[srcstr+2]*((0x10000-(grossx&0xffff))>>8)*
                                ((        (grossy&0xffff))>>8);
          sumX += ptr[srcstr+5]*((        (grossx&0xffff))>>8)*
                                ((        (grossy&0xffff))>>8);
          dst[0] = (unsigned char)((sumX/3)>>16);
        }
        else
          dst[0] = 0;
        dst++;
        grossx += xdX;
        grossy += ydX;
      }
      dst += add2dst;
    }
    return ERR_OK;
}

// YUV422 conversion
// rotate and resize ROI with bi-linear interpolation
MIA_RESULT_CODE IPL_GEOM_AffineTransformBLI_YUV422toGRAY(
	unsigned char* dst,         // destination image
	unsigned int dststr,        // stride in bytes
	const unsigned char* src,   // source image
	unsigned int srcstr,        // stride in bytes
	unsigned int xs,            // source image size
	unsigned int ys,            //
	double x,                   // source ROI
	double y,                   //
	const double* A,            // affine transformation matrix
	unsigned int W,             // destination ROI size
	unsigned int H)             //
{
	/*register*/ int grossx,grossy;
	int xdX,ydX,xdY,ydY;
	unsigned int sumX,i,j,add2dst;
	const unsigned char* ptr;

	// check arguments
	if ((src==NULL)||(dst==NULL)||(A==NULL))
		return ERR_GEN_NULLPOINTER;
	add2dst = dststr-W;
	// form deltas
	xdX = (int)(A[0]*0x00010000+.5);// (int)((cos(a)*w*0x00010000)/W+.5);
	ydX = (int)(A[2]*0x00010000+.5);// (int)((sin(a)*w*0x00010000)/W+.5);
	xdY = (int)(A[1]*0x00010000+.5);//-(int)((sin(a)*h*0x00010000)/H+.5);
	ydY = (int)(A[3]*0x00010000+.5);// (int)((cos(a)*h*0x00010000)/H+.5);
	// main cycle
	for (i=0;i<H;i++)
	{
		grossx = xdY*i+((int)(x*0x00010000));
		grossy = ydY*i+((int)(y*0x00010000));
		for (j=0;j<W;j++)
		{
			if ((grossx>=0)&&
				((grossx>>16)+1<(int)xs)&&
				(grossy>=0)&&
				((grossy>>16)+1<(int)ys) )
			{
				ptr = src+ ((grossy>>16)*srcstr)+(grossx>>16)*2;
				sumX = 0x8000;
				// Y
				sumX += ptr[0]*       ((0x10000-(grossx&0xffff))>>8)*
					((0x10000-(grossy&0xffff))>>8);
				sumX += ptr[2]*       ((        (grossx&0xffff))>>8)*
					((0x10000-(grossy&0xffff))>>8);
				sumX += ptr[srcstr+0]*((0x10000-(grossx&0xffff))>>8)*
					((        (grossy&0xffff))>>8);
				sumX += ptr[srcstr+2]*((        (grossx&0xffff))>>8)*
					((        (grossy&0xffff))>>8);

				dst[0] = (unsigned char)((sumX)>>16);
			}
			else
				dst[0] = 0;
			dst++;
			grossx += xdX;
			grossy += ydX;
		}
		dst += add2dst;
	}
	return ERR_OK;
}


// BAYER conversion
// R G R G R G
// G B G B G B
// rotate and resize ROI with bi-linear interpolation from Bayer RGGB
MIA_RESULT_CODE IPL_GEOM_AffineTransformBLI_BAYERtoGRAY(
	unsigned char* dst,         // destination image
	unsigned int dststr,        // stride in bytes
	const unsigned char* src,   // source image
	unsigned int srcstr,        // stride in bytes
	unsigned int xs,            // source image size
	unsigned int ys,            //
	double x,                   // source ROI
	double y,                   //
	const double* A,            // affine transformation matrix
	unsigned int W,             // destination ROI size
	unsigned int H)             //
{
/*register*/ int grossx,grossy;
  int xdX,ydX,xdY,ydY;
  unsigned int sumX,i,j,add2dst;
  const unsigned char* ptr;

    // check arguments
    if ((src==NULL)||(dst==NULL)||(A==NULL))
      return ERR_GEN_NULLPOINTER;
    add2dst = dststr-W;
    // form deltas
    xdX = (int)(A[0]/2*0x00010000+.5);// (int)((cos(a)*w*0x00010000)/W+.5);
    ydX = (int)(A[2]/2*0x00010000+.5);// (int)((sin(a)*w*0x00010000)/W+.5);
    xdY = (int)(A[1]/2*0x00010000+.5);//-(int)((sin(a)*h*0x00010000)/H+.5);
    ydY = (int)(A[3]/2*0x00010000+.5);// (int)((cos(a)*h*0x00010000)/H+.5);
    // main cycle
    for (i=0;i<H;i++)
    {
      grossx = xdY*i+((int)(x*0x00010000));
      grossy = ydY*i+((int)(y*0x00010000));
      for (j=0;j<W;j++)
      {
        if ((grossx>=0)&&
            ((grossx>>16)+1<(int)xs/2)&&
            (grossy>=0)&&
            ((grossy>>16)+1<(int)ys/2) )
        {
          ptr = src+ (((grossy>>16)*srcstr*2)+(grossx>>16)*2);
          sumX = 0x8000*3;
          // R
          sumX += ptr[0]*       ((0x10000-(grossx&0xffff))>>8)*
                                ((0x10000-(grossy&0xffff))>>8);
          sumX += ptr[2]*       ((        (grossx&0xffff))>>8)*
                                ((0x10000-(grossy&0xffff))>>8);
          sumX += ptr[2*srcstr+0]*((0x10000-(grossx&0xffff))>>8)*
                                ((        (grossy&0xffff))>>8);
          sumX += ptr[2*srcstr+2]*((        (grossx&0xffff))>>8)*
                                ((        (grossy&0xffff))>>8);
          // G
          sumX += ptr[1]*       ((0x10000-(grossx&0xffff))>>8)*
                                ((0x10000-(grossy&0xffff))>>8);
          sumX += ptr[3]*       ((        (grossx&0xffff))>>8)*
                                ((0x10000-(grossy&0xffff))>>8);
          sumX += ptr[2*srcstr+1]*((0x10000-(grossx&0xffff))>>8)*
                                ((        (grossy&0xffff))>>8);
          sumX += ptr[2*srcstr+3]*((        (grossx&0xffff))>>8)*
                                ((        (grossy&0xffff))>>8);
          // B
          sumX += ptr[srcstr+1]*       ((0x10000-(grossx&0xffff))>>8)*
                                ((0x10000-(grossy&0xffff))>>8);
          sumX += ptr[srcstr+3]*       ((        (grossx&0xffff))>>8)*
                                ((0x10000-(grossy&0xffff))>>8);
          sumX += ptr[3*srcstr+1]*((0x10000-(grossx&0xffff))>>8)*
                                ((        (grossy&0xffff))>>8);
          sumX += ptr[3*srcstr+3]*((        (grossx&0xffff))>>8)*
                                ((        (grossy&0xffff))>>8);
          dst[0] = (unsigned char)((sumX/3)>>16);
        }
        else
          dst[0] = 0;
        dst++;
        grossx += xdX;
        grossy += ydX;
      }
      dst += add2dst;
    }
    return ERR_OK;
}

// BAYER conversion
// R G R G R G
// G B G B G B
// rotate and resize ROI with nearest neighbour interpolation from Bayer RGGB
MIA_RESULT_CODE IPL_GEOM_AffineTransform_BAYERtoGRAYcrop(
	unsigned char* dst,         // destination image
	unsigned int dststr,        // stride in bytes
	const unsigned char* src,   // source image
	unsigned int srcstr,        // stride in bytes
	unsigned int xs,            // source image size
	unsigned int ys,            //
	double x,                   // source ROI
	double y,                   //
	const double* A,            // affine transformation matrix
	unsigned int W,             // destination ROI size
	unsigned int H)             //
{
/*register*/ int grossx,grossy;
  int xdX,ydX,xdY,ydY;
  unsigned int sumX,i,j,add2dst;
  const unsigned char* ptr;

    // check arguments
    if ((src==NULL)||(dst==NULL)||(A==NULL))
      return ERR_GEN_NULLPOINTER;
    add2dst = dststr-W;
    // form deltas
    xdX = (int)(A[0]/2*0x00010000+.5);// (int)((cos(a)*w*0x00010000)/W+.5);
    ydX = (int)(A[2]/2*0x00010000+.5);// (int)((sin(a)*w*0x00010000)/W+.5);
    xdY = (int)(A[1]/2*0x00010000+.5);//-(int)((sin(a)*h*0x00010000)/H+.5);
    ydY = (int)(A[3]/2*0x00010000+.5);// (int)((cos(a)*h*0x00010000)/H+.5);
    // main cycle
    for (i=0;i<H;i++)
    {
      grossx = xdY*i+((int)(x*0x00010000));
      grossy = ydY*i+((int)(y*0x00010000));
      for (j=0;j<W;j++)
      {
        if ((grossx>=0)&&
            ((grossx>>16)+1<(int)xs/2)&&
            (grossy>=0)&&
            ((grossy>>16)+1<(int)ys/2) )
        {
          ptr = src+ (((grossy>>16)*srcstr*2)+(grossx>>16)*2);
          sumX = 0x8000*3;
          // R
          sumX += ptr[0];
          // G
          sumX += ptr[1];
          // B
          sumX += ptr[srcstr+1];
          dst[0] = (unsigned char)((sumX/3));
        }
        else
          dst[0] = 0;
        dst++;
        grossx += xdX;
        grossy += ydX;
      }
      dst += add2dst;
    }
    return ERR_OK;
}

// rotate and resize ROI with nearest-neighbor interpolation
MIA_RESULT_CODE IPL_GEOM_AffineTransformNNI_RGB(
  unsigned char* dst,         // destination image
  unsigned int dststr,        // stride in bytes
  const unsigned char* src,   // source image
  unsigned int srcstr,        // stride in bytes
  unsigned int xs,            // source image size
  unsigned int ys,            //
  double x,                   // source ROI
  double y,                   //
  const double* A,            // affine transformation matrix
  unsigned int W,             // destination ROI size
  unsigned int H)             //
{
  /*register*/ int grossx,grossy;
  int xdX,ydX,xdY,ydY;
  unsigned int i,j,add2dst;
  const unsigned char* ptr;

    // check arguments
    if ((src==NULL)||(dst==NULL)||(A==NULL))
      return ERR_GEN_NULLPOINTER;
    add2dst = dststr-W*3;
    // form deltas
    xdX = (int)(A[0]*0x00010000+.5);// (int)((cos(a)*w*0x00010000)/W+.5);
    ydX = (int)(A[2]*0x00010000+.5);// (int)((sin(a)*w*0x00010000)/W+.5);
    xdY = (int)(A[1]*0x00010000+.5);//-(int)((sin(a)*h*0x00010000)/H+.5);
    ydY = (int)(A[3]*0x00010000+.5);// (int)((cos(a)*h*0x00010000)/H+.5);
    // main cycle
    for (i=0;i<H;i++)
    {
      grossx = xdY*i+((int)(x*0x00010000));
      grossy = ydY*i+((int)(y*0x00010000));
      for (j=0;j<W;j++)
      {
        if ((grossx>=0)&&
            ((grossx>>16)+1<(int)xs)&&
            (grossy>=0)&&
            ((grossy>>16)+1<(int)ys) )
        {
          ptr = src+ ((grossy>>16)*srcstr)+(grossx>>16)*3;
          if ((grossx&0xffff)>0x8000)
            ptr += 3;
          if ((grossy&0xffff)>0x8000)
            ptr += srcstr;
          dst[0] = ptr[0];
          dst[1] = ptr[1];
          dst[2] = ptr[2];
        }
        else
          dst[0] = dst[1] = dst[2] = 0;
        dst += 3;
        grossx += xdX;
        grossy += ydX;
      }
      dst += add2dst;
    }
    return ERR_OK;
}

// mirrors image relative to y axis (in-place transformation)
void IPL_GEOM_MirrorX3(
  unsigned char* im,        // destination
  unsigned int xs,         // image size
  unsigned int ys)         //
{
  unsigned int *ptr1,*ptr2,val_s,val_d;

    ptr1 = (unsigned int*)im;
    xs /= 8;
    while (ys--)
    {
      ptr2 = ptr1+xs*2-1;
      do
      {
        // beginning
        val_d = val_s = *ptr1;
        val_s >>= 8;
        val_d <<= 8;
        val_d += (val_s&0xff);
        val_s >>= 8;
        val_d <<= 8;
        val_d += (val_s&0xff);
        val_d <<= 8;
        val_d += (val_s>>8);
        // this can be easily done by XCHG !!!
        val_s = *ptr2;
        *ptr2 = val_d;
        val_d = val_s;
        // ending
        val_s >>= 8;
        val_d <<= 8;
        val_d += (val_s&0xff);
        val_s >>= 8;
        val_d <<= 8;
        val_d += (val_s&0xff);
        val_d <<= 8;
        val_d += (val_s>>8);
        *ptr1 = val_d;
        ptr1++;
        ptr2--;
      }
      while (ptr2>ptr1);
      ptr1 += xs;
    }
}

/* used in DSP face identification earlier
// mirrors image relative to y axis (in-place transformation)
void IPL_GEOM_MirrorX3(
  unsigned char* im,        // destination
  unsigned int xs,         // image size
  unsigned int ys)         //
{
  unsigned char val;
  unsigned char* km;

    km = im+xs-1;
    while (ys--)
    {
      while ((unsigned int)im<(unsigned int)km)
      {
        val = *im;
        *im = *km;
        *km = val;
        im++;
        km--;
      }
      im += (xs+1)/2;
      km +=  xs+xs/2;
    }
}
*/

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

// rotate and resize ROI with bi-linear interpolation
MIA_RESULT_CODE IPL_GEOM_AffineTransformBLI(
  unsigned char* dst,         // destination image
  unsigned int dststr,        // stride in bytes
  const unsigned char* src,   // source image
  unsigned int srcstr,        // stride in bytes
  unsigned int xs,            // source image size
  unsigned int ys,            //
  double x,                   // corner of source ROI
  double y,                   //
  const double* A,            // affine transformation matrix
  unsigned int W,             // destination ROI size
  unsigned int H)             //
{

  /*register*/ int grossx,grossy;
  int xdX,ydX,xdY,ydY;
  unsigned int sum,i,j,add2dst;
  const unsigned char* ptr;

    // check arguments
    if ((src==NULL)||(dst==NULL)||(A==NULL))
      return ERR_GEN_NULLPOINTER;
    // +++ here one should insert the separating and different treatment
    //     of three cases:
    //      1) no rotation, no resize - simple copying ROI
    //      2) no rotation - simple resize
    //      3) rotation is present - this is the most universal most complex
    //         and thus slowest method. It is implemented now.
    //         Furthermore, rotations by degrees of pi/2 without resising
    //         can be outgrouped to one more case - mirroring.
    // ---
    add2dst = dststr-W;
    // form deltas
    xdX = (int)(A[0]*0x00010000+.5);// (int)((cos(a)*w*0x00010000)/W+.5);
    ydX = (int)(A[2]*0x00010000+.5);// (int)((sin(a)*w*0x00010000)/W+.5);
    xdY = (int)(A[1]*0x00010000+.5);//-(int)((sin(a)*h*0x00010000)/H+.5);
    ydY = (int)(A[3]*0x00010000+.5);// (int)((cos(a)*h*0x00010000)/H+.5);
    // main cycle
    for (i=0;i<H;i++)
    {
      grossx = xdY*i+((int)(x*0x00010000));
      grossy = ydY*i+((int)(y*0x00010000));
      for (j=0;j<W;j++)
      {
        if ((grossx>=0)&&
            ((grossx>>16)+1<(int)xs)&&
            (grossy>=0)&&
            ((grossy>>16)+1<(int)ys) )
        {
          sum = 0x8000;
          ptr = src+ ((grossy>>16)*srcstr)+(grossx>>16);
          sum += ptr[0]*        ((0x10000-(grossx&0xffff))>>8)*
                                ((0x10000-(grossy&0xffff))>>8);
          sum += ptr[1]*        ((        (grossx&0xffff))>>8)*
                                ((0x10000-(grossy&0xffff))>>8);
          sum += ptr[srcstr]*   ((0x10000-(grossx&0xffff))>>8)*
                                ((        (grossy&0xffff))>>8);
          sum += ptr[srcstr+1]* ((        (grossx&0xffff))>>8)*
                                ((        (grossy&0xffff))>>8);
          *dst = (unsigned char)(sum>>16);
        }
        else
          *dst = 0;
        dst++;
        grossx += xdX;
        grossy += ydX;
      }
      dst += add2dst;
    }
    return ERR_OK;
}

/*
// reference function for RotateResizeBLI - straightforward realization
MIA_RESULT_CODE IPL_RotateResizeBLI_Reference(
          unsigned char* dst,     // destination image
          unsigned int dststr,   // stride in bytes
    const unsigned char* src,     // source image
          unsigned int srcstr,   // stride in bytes
          unsigned int xs,       // source image size
          unsigned int ys,       //
          double x,               // source ROI
          double y,               //
          double w,               //
          double h,               //
          double a,               // angle of sorce ROI tilt
          unsigned int W,        // destination ROI size
          unsigned int H)        //
{
  unsigned int i,j;
  double yd,xd,yp,xp,wsin,wcos,hsin,hcos;
  int xi,yi;

    // check arguments
    if ((src==NULL)||(dst==NULL))
      return ERR_GEN_NULLPOINTER;
    wsin = w*sin(a)/W;
    hcos = h*cos(a)/H;
    wcos = w*cos(a)/W;
    hsin = h*sin(a)/H;
    for (i=0;i<H;i++)
    {
      for (j=0;j<W;j++)
      {
        yd = y+(j*wsin)+(i*hcos);
        yi = int(yd);
        yp = yd-yi;
        xd = x+(j*wcos)-(i*hsin);
        xi = int(xd);
        xp = xd-xi;
        if ((yi>0)&&(yi<(int)ys-1)&&(xi>0)&&(xi<(int)xs-1))
          dst[i*dststr+j] = (unsigned char)(
              src[ yi   *srcstr+ xi   ]*(1.-xp)*(1.-yp)+
              src[ yi   *srcstr+(xi+1)]*    xp *(1.-yp)+
              src[(yi+1)*srcstr+ xi   ]*(1.-xp)*    yp +
              src[(yi+1)*srcstr+(xi+1)]*    xp *    yp );
        else
          dst[i*dststr+j] = 0;
      }
    }
    return ERR_OK;
}
*/

// mirror image relative to vertical axis
MIA_RESULT_CODE IPL_GEOM_MirrorX(
  unsigned char* dst,           // destination
  unsigned int dststr,          // stride
  const unsigned char* src,     // source
  unsigned int srcstr,          // stride
  unsigned int w,               // ROI size
  unsigned int h)               //
{
  unsigned int i,j;

    // check arguments
    if ((dst==NULL)||(src==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((dststr==0)||(srcstr==0)||(w==0)||(h==0))
      return ERR_GEN_NOMEMORY;
    for (i=0;i<h;i++)
      for (j=0;j<w;j++)
        dst[i*dststr+j] = src[i*srcstr+w-j-1];
    return ERR_OK;
}

// mirror image relative to vertical axis
// different BpP values possible
MIA_RESULT_CODE IPL_GEOM_MirrorX2(
  unsigned char* dst,           // destination
  const unsigned char* src,     // source
  int xs,                       // image size
  int ys,                       //
  int BpP)                      // bytes per pixel
{
  int i,j,k;

    for (i=0;i<ys;i++)
      for (j=0;j<xs;j++)
        for (k=0;k<BpP;k++)
          dst[((i)*xs+(j))*BpP+k] = src[((i)*xs+(xs-j-1))*BpP+k];
    return ERR_OK;
}

// mirror image relative to horizontal axis
// different BpP values possible
MIA_RESULT_CODE IPL_GEOM_MirrorY(
  unsigned char* dst,           // destination
  const unsigned char* src,     // source
  int xs,                       // image size
  int ys,                       //
  int BpP)                      // bytes per pixel
{
  int i,j,k;

    for (j=0;j<xs;j++)
      for (i=0;i<ys;i++)
        for (k=0;k<BpP;k++)
          dst[((i)*xs+(j))*BpP+k] = src[((ys-i-1)*xs+(j))*BpP+k];
    return ERR_OK;
}

// flip image 90 degrees counterclockwise
// BpP!=1 can be used
MIA_RESULT_CODE IPL_GEOM_Flip90(
  unsigned char* dst,           // destination
  const unsigned char* src,     // source
  int xs,                       // image size
  int ys,                       //
  int BpP)                      // bytes per pixel
{
  int i,j,k;

    for (j=0;j<xs;j++)
      for (i=0;i<ys;i++)
        for (k=0;k<BpP;k++)
          dst[((xs-j-1)*ys+(i))*BpP+k] = src[((i)*xs+(j))*BpP+k];
    return ERR_OK;
}

// using non-concentric pupil and iris circles
MIA_RESULT_CODE IPL_GEOM_GetPolarImage1(  // MIA 11.08.04
  unsigned char* dst, // OUT: polar data block 
  int phi,            // IN: number of lines in polar image
  int rho,            // IN: number of columns in polar image
  const unsigned char* src,// IN:  eye image
  int xs,             // IN: image size
  int ys,             //
  int xp,             // center coordinates of pupil (inner boundary)
  int yp,             //
  int rp,             // radius of pupil
  int xi,             // center coordinates of iris (outer boundary)
  int yi,             //
  int ri)             // radius of iris
{
  int rhoidx,phiidx,xint,yint;
  double f,alpha,xbeg,ybeg,xend,yend,xstep,ystep,tmpdbl,xfrac,yfrac,val;

    for (phiidx=0;phiidx<phi;phiidx++)
    {
      alpha = phiidx*2*PI/phi;
      xbeg = xp+rp*cos(alpha);
      ybeg = yp+rp*sin(alpha);
      f = (xp-xi)*sin(alpha)-(yp-yi)*cos(alpha);
      tmpdbl = ri*ri-f;
      if (tmpdbl<=0.)
        return ERR_GEN_INTERNAL;
      tmpdbl = sqrt(tmpdbl);
      xend = xi+f*sin(alpha)+tmpdbl*cos(alpha);
      yend = yi+f*cos(alpha)+tmpdbl*sin(alpha);
      xstep = (xend-xbeg)/rho;
      ystep = (yend-ybeg)/rho;
      xbeg -= xstep/2.;
      ybeg -= ystep/2.;
      for (rhoidx=0;rhoidx<rho;rhoidx++)
      {
        xbeg += xstep;
        ybeg += ystep;
        if ((xbeg<0)||(xbeg>=xs-1)||(ybeg<0)||(ybeg>=ys-1))
        {
          dst[(rho-rhoidx-1)*phi+phiidx] = 0;
          continue;
        }
        xint = (int)xbeg;
        xfrac = xbeg-xint;
        yint = (int)ybeg;
        yfrac = ybeg-yint;
        val = (1.-xfrac)*(1.-yfrac)*src[ yint   *xs+ xint   ]+
                  xfrac *(1.-yfrac)*src[ yint   *xs+(xint+1)]+
              (1.-xfrac)*    yfrac *src[(yint+1)*xs+ xint   ]+
                  xfrac *    yfrac *src[(yint+1)*xs+(xint+1)];
        if (val<0)
          val = 0;
        if (val>255)
          val = 255;
        dst[(rho-rhoidx-1)*phi+phiidx] = (unsigned char)val;
      }
    }
    return ERR_OK;
}
