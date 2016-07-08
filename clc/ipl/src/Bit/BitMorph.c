/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  26 November 2001                                       */
/*      Revision: 1.0.00                                                 */
/*      Purpose:  bit image operations                                   */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#define __FILENUM__ 3 // __FILENUM__TAG3

#include <memory.h>
#include "ipl.h"

// morphological erosion with circle 9x9 window
staticFunc void BitMorph_Erode4(
            unsigned char* dst,     // destination
      const unsigned char* src,     // source
            unsigned int w,        // width in pixels, aligned by DWORD
            unsigned int h)        // height
{
  unsigned i,j,ww;
  unsigned int* ptr;

    w /= 8;   // w now is in bytes
    ww = w/4;
    memset(dst,0,w*h);
    for (i=4;i<h-4;i++)
      for (j=0;j<w;j+=3)
      {
        ptr = (unsigned int*)&src[ i   *w+j];
        *((unsigned int*)&dst[ i   *w+j]) |=
                (                                                 (ptr[-4*ww]<<1)&(ptr[-4*ww]   )&(ptr[-4*ww]>>1)&
                                  (ptr[-3*ww]<<3)&(ptr[-3*ww]<<2)&(ptr[-3*ww]<<1)&(ptr[-3*ww]   )&(ptr[-3*ww]>>1)&(ptr[-3*ww]>>2)&(ptr[-3*ww]>>3)&
                                  (ptr[-2*ww]<<3)&(ptr[-2*ww]<<2)&(ptr[-2*ww]<<1)&(ptr[-2*ww]   )&(ptr[-2*ww]>>1)&(ptr[-2*ww]>>2)&(ptr[-2*ww]>>3)&
                  (ptr[-1*ww]<<4)&(ptr[-1*ww]<<3)&(ptr[-1*ww]<<2)&(ptr[-1*ww]<<1)&(ptr[-1*ww]   )&(ptr[-1*ww]>>1)&(ptr[-1*ww]>>2)&(ptr[-1*ww]>>3)&(ptr[-1*ww]>>4)&
                  (ptr[    0]<<4)&(ptr[    0]<<3)&(ptr[    0]<<2)&(ptr[    0]<<1)&(ptr[    0]   )&(ptr[    0]>>1)&(ptr[    0]>>2)&(ptr[    0]>>3)&(ptr[    0]>>4)&
                  (ptr[ 1*ww]<<4)&(ptr[ 1*ww]<<3)&(ptr[ 1*ww]<<2)&(ptr[ 1*ww]<<1)&(ptr[ 1*ww]   )&(ptr[ 1*ww]>>1)&(ptr[ 1*ww]>>2)&(ptr[ 1*ww]>>3)&(ptr[ 1*ww]>>4)&
                                  (ptr[ 2*ww]<<3)&(ptr[ 2*ww]<<2)&(ptr[ 2*ww]<<1)&(ptr[ 2*ww]   )&(ptr[ 2*ww]>>1)&(ptr[ 2*ww]>>2)&(ptr[ 2*ww]>>3)&
                                  (ptr[ 3*ww]<<3)&(ptr[ 3*ww]<<2)&(ptr[ 3*ww]<<1)&(ptr[ 3*ww]   )&(ptr[ 3*ww]>>1)&(ptr[ 3*ww]>>2)&(ptr[ 3*ww]>>3)&
                                                                  (ptr[ 4*ww]<<1)&(ptr[ 4*ww]   )&(ptr[ 4*ww]>>1)                                                  )
                & 0x0ffffff0;
      }
    if ( (w-4)%3 )
      for (i=4;i<h-4;i++)
      {
        j = 0; // ?
        ptr = (unsigned int*)&src[ i   *w+j];
        *((unsigned int*)&dst[(i+1)*w-4]) |=
                (                                                 (ptr[-4*ww]<<1)&(ptr[-4*ww]   )&(ptr[-4*ww]>>1)&
                                  (ptr[-3*ww]<<3)&(ptr[-3*ww]<<2)&(ptr[-3*ww]<<1)&(ptr[-3*ww]   )&(ptr[-3*ww]>>1)&(ptr[-3*ww]>>2)&(ptr[-3*ww]>>3)&
                                  (ptr[-2*ww]<<3)&(ptr[-2*ww]<<2)&(ptr[-2*ww]<<1)&(ptr[-2*ww]   )&(ptr[-2*ww]>>1)&(ptr[-2*ww]>>2)&(ptr[-2*ww]>>3)&
                  (ptr[-1*ww]<<4)&(ptr[-1*ww]<<3)&(ptr[-1*ww]<<2)&(ptr[-1*ww]<<1)&(ptr[-1*ww]   )&(ptr[-1*ww]>>1)&(ptr[-1*ww]>>2)&(ptr[-1*ww]>>3)&(ptr[-1*ww]>>4)&
                  (ptr[    0]<<4)&(ptr[    0]<<3)&(ptr[    0]<<2)&(ptr[    0]<<1)&(ptr[    0]   )&(ptr[    0]>>1)&(ptr[    0]>>2)&(ptr[    0]>>3)&(ptr[    0]>>4)&
                  (ptr[ 1*ww]<<4)&(ptr[ 1*ww]<<3)&(ptr[ 1*ww]<<2)&(ptr[ 1*ww]<<1)&(ptr[ 1*ww]   )&(ptr[ 1*ww]>>1)&(ptr[ 1*ww]>>2)&(ptr[ 1*ww]>>3)&(ptr[ 1*ww]>>4)&
                                  (ptr[ 2*ww]<<3)&(ptr[ 2*ww]<<2)&(ptr[ 2*ww]<<1)&(ptr[ 2*ww]   )&(ptr[ 2*ww]>>1)&(ptr[ 2*ww]>>2)&(ptr[ 2*ww]>>3)&
                                  (ptr[ 3*ww]<<3)&(ptr[ 3*ww]<<2)&(ptr[ 3*ww]<<1)&(ptr[ 3*ww]   )&(ptr[ 3*ww]>>1)&(ptr[ 3*ww]>>2)&(ptr[ 3*ww]>>3)&
                                                                  (ptr[ 4*ww]<<1)&(ptr[ 4*ww]   )&(ptr[ 4*ww]>>1)                                                  )
                & 0x0ffffff0;
      }
}

// morphological erosion with circle 7x7 window
staticFunc void BitMorph_Erode3(
            unsigned char* dst,     // destination
      const unsigned char* src,     // source
            unsigned int w,        // width in pixels, aligned by DWORD
            unsigned int h)        // height
{
  unsigned i,j,ww;
  unsigned int* ptr;

    w /= 8;   // w now is in bytes
    ww = w/4;
    memset(dst,0,w*h);
    for (i=3;i<h-3;i++)
      for (j=0;j<w;j+=3)
      {
        ptr = (unsigned int*)&src[ i   *w+j];
        *((unsigned int*)&dst[ i   *w+j]) |=
                (                                 (ptr[-3*ww]<<1)&(ptr[-3*ww]   )&(ptr[-3*ww]>>1)&
                                  (ptr[-2*ww]<<2)&(ptr[-2*ww]<<1)&(ptr[-2*ww]   )&(ptr[-2*ww]>>1)&(ptr[-2*ww]>>2)&
                  (ptr[-1*ww]<<3)&(ptr[-1*ww]<<2)&(ptr[-1*ww]<<1)&(ptr[-1*ww]   )&(ptr[-1*ww]>>1)&(ptr[-1*ww]>>2)&(ptr[-1*ww]>>3)&
                  (ptr[    0]<<3)&(ptr[    0]<<2)&(ptr[    0]<<1)&(ptr[    0]   )&(ptr[    0]>>1)&(ptr[    0]>>2)&(ptr[    0]>>3)&
                  (ptr[ 1*ww]<<3)&(ptr[ 1*ww]<<2)&(ptr[ 1*ww]<<1)&(ptr[ 1*ww]   )&(ptr[ 1*ww]>>1)&(ptr[ 1*ww]>>2)&(ptr[ 1*ww]>>3)&
                                  (ptr[ 2*ww]<<2)&(ptr[ 2*ww]<<1)&(ptr[ 2*ww]   )&(ptr[ 2*ww]>>1)&(ptr[ 2*ww]>>2)&
                                                  (ptr[ 3*ww]<<1)&(ptr[ 3*ww]   )&(ptr[ 3*ww]>>1)                                 )
                & 0x0ffffff0;
      }
    if ( (w-4)%3 )
      for (i=3;i<h-3;i++)
      {
        j = 0; // ?
        ptr = (unsigned int*)&src[ i   *w+j];
        *((unsigned int*)&dst[(i+1)*w-4]) |=
                (                                 (ptr[-3*ww]<<1)&(ptr[-3*ww]   )&(ptr[-3*ww]>>1)&
                                  (ptr[-2*ww]<<2)&(ptr[-2*ww]<<1)&(ptr[-2*ww]   )&(ptr[-2*ww]>>1)&(ptr[-2*ww]>>2)&
                  (ptr[-1*ww]<<3)&(ptr[-1*ww]<<2)&(ptr[-1*ww]<<1)&(ptr[-1*ww]   )&(ptr[-1*ww]>>1)&(ptr[-1*ww]>>2)&(ptr[-1*ww]>>3)&
                  (ptr[    0]<<3)&(ptr[    0]<<2)&(ptr[    0]<<1)&(ptr[    0]   )&(ptr[    0]>>1)&(ptr[    0]>>2)&(ptr[    0]>>3)&
                  (ptr[ 1*ww]<<3)&(ptr[ 1*ww]<<2)&(ptr[ 1*ww]<<1)&(ptr[ 1*ww]   )&(ptr[ 1*ww]>>1)&(ptr[ 1*ww]>>2)&(ptr[ 1*ww]>>3)&
                                  (ptr[ 2*ww]<<2)&(ptr[ 2*ww]<<1)&(ptr[ 2*ww]   )&(ptr[ 2*ww]>>1)&(ptr[ 2*ww]>>2)&
                                                  (ptr[ 3*ww]<<1)&(ptr[ 3*ww]   )&(ptr[ 3*ww]>>1)                                 )
                & 0x0ffffff0;
      }
}

// morphological dilation with circle 9x9 window
staticFunc void BitMorph_Dilate4(
            unsigned char* dst,     // destination
      const unsigned char* src,     // source
            unsigned int w,        // width in pixels, aligned by DWORD
            unsigned int h)        // height
{
  unsigned i,j,ww;
  unsigned int* ptr;

    w /= 8;     // w now is in bytes
    ww = w/4;
    memset(dst,0,w*h);
    for (i=4;i<h-4;i++)
      for (j=0;j<w;j+=3)
      {
        ptr = (unsigned int*)&src[ i   *w+j];
        *((unsigned int*)&dst[ i   *w+j]) |=
                (                                                 (ptr[-4*ww]<<1)|(ptr[-4*ww]   )|(ptr[-4*ww]>>1)|
                                  (ptr[-3*ww]<<3)|(ptr[-3*ww]<<2)|(ptr[-3*ww]<<1)|(ptr[-3*ww]   )|(ptr[-3*ww]>>1)|(ptr[-3*ww]>>2)|(ptr[-3*ww]>>3)|
                                  (ptr[-2*ww]<<3)|(ptr[-2*ww]<<2)|(ptr[-2*ww]<<1)|(ptr[-2*ww]   )|(ptr[-2*ww]>>1)|(ptr[-2*ww]>>2)|(ptr[-2*ww]>>3)|
                  (ptr[-1*ww]<<4)|(ptr[-1*ww]<<3)|(ptr[-1*ww]<<2)|(ptr[-1*ww]<<1)|(ptr[-1*ww]   )|(ptr[-1*ww]>>1)|(ptr[-1*ww]>>2)|(ptr[-1*ww]>>3)|(ptr[-1*ww]>>4)|
                  (ptr[    0]<<4)|(ptr[    0]<<3)|(ptr[    0]<<2)|(ptr[    0]<<1)|(ptr[    0]   )|(ptr[    0]>>1)|(ptr[    0]>>2)|(ptr[    0]>>3)|(ptr[    0]>>4)|
                  (ptr[ 1*ww]<<4)|(ptr[ 1*ww]<<3)|(ptr[ 1*ww]<<2)|(ptr[ 1*ww]<<1)|(ptr[ 1*ww]   )|(ptr[ 1*ww]>>1)|(ptr[ 1*ww]>>2)|(ptr[ 1*ww]>>3)|(ptr[ 1*ww]>>4)|
                                  (ptr[ 2*ww]<<3)|(ptr[ 2*ww]<<2)|(ptr[ 2*ww]<<1)|(ptr[ 2*ww]   )|(ptr[ 2*ww]>>1)|(ptr[ 2*ww]>>2)|(ptr[ 2*ww]>>3)|
                                  (ptr[ 3*ww]<<3)|(ptr[ 3*ww]<<2)|(ptr[ 3*ww]<<1)|(ptr[ 3*ww]   )|(ptr[ 3*ww]>>1)|(ptr[ 3*ww]>>2)|(ptr[ 3*ww]>>3)|
                                                                  (ptr[ 4*ww]<<1)|(ptr[ 4*ww]   )|(ptr[ 4*ww]>>1)                                                 )
                & 0x0ffffff0;
      }
    if ( (w-4)%3 )
      for (i=4;i<h-4;i++)
      {
        j = 0; // ?
        ptr = (unsigned int*)&src[ i   *w+j];
        *((unsigned int*)&dst[(i+1)*w-4]) |=
                (                                                 (ptr[-4*ww]<<1)|(ptr[-4*ww]   )|(ptr[-4*ww]>>1)|
                                  (ptr[-3*ww]<<3)|(ptr[-3*ww]<<2)|(ptr[-3*ww]<<1)|(ptr[-3*ww]   )|(ptr[-3*ww]>>1)|(ptr[-3*ww]>>2)|(ptr[-3*ww]>>3)|
                                  (ptr[-2*ww]<<3)|(ptr[-2*ww]<<2)|(ptr[-2*ww]<<1)|(ptr[-2*ww]   )|(ptr[-2*ww]>>1)|(ptr[-2*ww]>>2)|(ptr[-2*ww]>>3)|
                  (ptr[-1*ww]<<4)|(ptr[-1*ww]<<3)|(ptr[-1*ww]<<2)|(ptr[-1*ww]<<1)|(ptr[-1*ww]   )|(ptr[-1*ww]>>1)|(ptr[-1*ww]>>2)|(ptr[-1*ww]>>3)|(ptr[-1*ww]>>4)|
                  (ptr[    0]<<4)|(ptr[    0]<<3)|(ptr[    0]<<2)|(ptr[    0]<<1)|(ptr[    0]   )|(ptr[    0]>>1)|(ptr[    0]>>2)|(ptr[    0]>>3)|(ptr[    0]>>4)|
                  (ptr[ 1*ww]<<4)|(ptr[ 1*ww]<<3)|(ptr[ 1*ww]<<2)|(ptr[ 1*ww]<<1)|(ptr[ 1*ww]   )|(ptr[ 1*ww]>>1)|(ptr[ 1*ww]>>2)|(ptr[ 1*ww]>>3)|(ptr[ 1*ww]>>4)|
                                  (ptr[ 2*ww]<<3)|(ptr[ 2*ww]<<2)|(ptr[ 2*ww]<<1)|(ptr[ 2*ww]   )|(ptr[ 2*ww]>>1)|(ptr[ 2*ww]>>2)|(ptr[ 2*ww]>>3)|
                                  (ptr[ 3*ww]<<3)|(ptr[ 3*ww]<<2)|(ptr[ 3*ww]<<1)|(ptr[ 3*ww]   )|(ptr[ 3*ww]>>1)|(ptr[ 3*ww]>>2)|(ptr[ 3*ww]>>3)|
                                                                  (ptr[ 4*ww]<<1)|(ptr[ 4*ww]   )|(ptr[ 4*ww]>>1)                                                 )
                & 0x0ffffff0;
      }
}

// morphological dilation with circle 7x7 window
staticFunc void BitMorph_Dilate3(
            unsigned char* dst,     // destination
      const unsigned char* src,     // source
            unsigned int w,        // width in pixels, aligned by DWORD
            unsigned int h)        // height
{
  unsigned i,j,ww;
  unsigned int* ptr;

    w /= 8;     // w now is in bytes
    ww = w/4;
    memset(dst,0,w*h);
    for (i=3;i<h-3;i++)
      for (j=0;j<w;j+=3)
      {
        ptr = (unsigned int*)&src[ i   *w+j];
        *((unsigned int*)&dst[ i   *w+j]) |=
                (                                 (ptr[-3*ww]<<1)|(ptr[-3*ww]   )|(ptr[-3*ww]>>1)|
                                  (ptr[-2*ww]<<2)|(ptr[-2*ww]<<1)|(ptr[-2*ww]   )|(ptr[-2*ww]>>1)|(ptr[-2*ww]>>2)|
                  (ptr[-1*ww]<<3)|(ptr[-1*ww]<<2)|(ptr[-1*ww]<<1)|(ptr[-1*ww]   )|(ptr[-1*ww]>>1)|(ptr[-1*ww]>>2)|(ptr[-1*ww]>>3)|
                  (ptr[    0]<<3)|(ptr[    0]<<2)|(ptr[    0]<<1)|(ptr[    0]   )|(ptr[    0]>>1)|(ptr[    0]>>2)|(ptr[    0]>>3)|
                  (ptr[ 1*ww]<<3)|(ptr[ 1*ww]<<2)|(ptr[ 1*ww]<<1)|(ptr[ 1*ww]   )|(ptr[ 1*ww]>>1)|(ptr[ 1*ww]>>2)|(ptr[ 1*ww]>>3)|
                                  (ptr[ 2*ww]<<2)|(ptr[ 2*ww]<<1)|(ptr[ 2*ww]   )|(ptr[ 2*ww]>>1)|(ptr[ 2*ww]>>2)|
                                                  (ptr[ 3*ww]<<1)|(ptr[ 3*ww]   )|(ptr[ 3*ww]>>1)                                 )
                & 0x0ffffff0;
      }
    if ( (w-4)%3 )
      for (i=3;i<h-3;i++)
      {
        j = 0; // ?
        ptr = (unsigned int*)&src[ i   *w+j];
        *((unsigned int*)&dst[(i+1)*w-4]) |=
                (                                 (ptr[-3*ww]<<1)|(ptr[-3*ww]   )|(ptr[-3*ww]>>1)|
                                  (ptr[-2*ww]<<2)|(ptr[-2*ww]<<1)|(ptr[-2*ww]   )|(ptr[-2*ww]>>1)|(ptr[-2*ww]>>2)|
                  (ptr[-1*ww]<<3)|(ptr[-1*ww]<<2)|(ptr[-1*ww]<<1)|(ptr[-1*ww]   )|(ptr[-1*ww]>>1)|(ptr[-1*ww]>>2)|(ptr[-1*ww]>>3)|
                  (ptr[    0]<<3)|(ptr[    0]<<2)|(ptr[    0]<<1)|(ptr[    0]   )|(ptr[    0]>>1)|(ptr[    0]>>2)|(ptr[    0]>>3)|
                  (ptr[ 1*ww]<<3)|(ptr[ 1*ww]<<2)|(ptr[ 1*ww]<<1)|(ptr[ 1*ww]   )|(ptr[ 1*ww]>>1)|(ptr[ 1*ww]>>2)|(ptr[ 1*ww]>>3)|
                                  (ptr[ 2*ww]<<2)|(ptr[ 2*ww]<<1)|(ptr[ 2*ww]   )|(ptr[ 2*ww]>>1)|(ptr[ 2*ww]>>2)|
                                                  (ptr[ 3*ww]<<1)|(ptr[ 3*ww]   )|(ptr[ 3*ww]>>1)                                 )
                & 0x0ffffff0;
      }
}

// bit image morphology
// WARNING: checking boundaries is on behalf of a caller !
MIA_RESULT_CODE IPL_BIT_Morphology(
  unsigned char* dst,             // destination
  unsigned int dststr,           // stride in bytes
  const unsigned char* src,       // source
  unsigned int srcstr,           // stride in bytes
  unsigned int xs,               // image size in pixels
  unsigned int ys,               //
  IPL_BIT_OPERATION operation,    // morphological operation
  unsigned int size)             // size of operation
{
    // check arguments
    if ((dst==NULL)||(src==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((dststr==0)||(srcstr==0)||(xs==0)||(ys==0))
      return ERR_GEN_NO_DATA;
    if ((xs!=srcstr*8)||(xs!=dststr*8))
      return ERR_GEN_NOTIMPL;
    if (xs&31)
      return ERR_GEN_INVALID_PARAMETER;
    //
    if (operation==IPL_BIT_ERODE)
    {
      if (size==3)
        BitMorph_Erode3(dst,src,xs,ys);
      if (size==4)
        BitMorph_Erode4(dst,src,xs,ys);
    }
    if (operation==IPL_BIT_DILATE)
    {
      if (size==3)
        BitMorph_Dilate3(dst,src,xs,ys);
      if (size==4)
        BitMorph_Dilate4(dst,src,xs,ys);
    }
    return ERR_OK;
}

MIA_RESULT_CODE IPL_BIT_ByteImageToBitImage(
  unsigned int* dst,        // dest image
  int dst_str,              // dest stride in DWORDs
  const unsigned char* src, // source image
  int src_str,              // source stride in bytes
  int width,                // width of source region to be converted
  int height)               // height
{
  int i,j;

    // check arguments
    if ((dst==NULL)||(src==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((dst_str*4<width)||(src_str<width))
      return ERR_GEN_INSUFFICIENT_BUFFER;
    if ((dst_str==0)||(src_str==0)||(height==0))
      return ERR_GEN_NO_DATA;
    //
    memset(dst,0,dst_str*height*4);
    // should be optimized...
    for (i=0;i<height;i++)
      for (j=0;j<width;j++)
        dst[i*dst_str+j/32] |= ((src[i*src_str+j])?(1<<(j%32)):0);
    return ERR_OK;
}
