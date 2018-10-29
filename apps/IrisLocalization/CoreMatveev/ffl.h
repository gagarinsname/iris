/*------------------------------------------------------------------------*/
/*                                                                        */
/*      Created:  14 March 2001 - primarily as video input library        */
/*      Modified: 14 April 2001 - file sequences added (version=2)        */
/*      Modified: 12 October 2001 - DAT files as sequences (version=3)    */
/*      Modified: 3 September 2002 - DTC file handling added (version=4)  */
/*      Modified: 5 March 2004 - Video input removed (version=5)          */
/*      Modified: 19 July 2010 - revised for CodCon (version=6)           */
/*      Modified: 6 August 2012 - DTC4 added                              */
/*      Revision: 6.0.00                                                  */
/*      Purpose:  File Format Library                                     */
/*      Authors:                                                          */
/*        Ivan Matveev                                                    */
/*                                                                        */
/*------------------------------------------------------------------------*/

#ifndef __ffl_h__
#define __ffl_h__

#include <stdio.h>
#include "errorcodes.h"
#include "stddefs.h"
#ifdef linux
  #include <sys/time.h>
#else
  #include <time.h>
#endif //linux

EXTERNC MIA_RESULT_CODE FFL_GetVersionInfo(
  char* pcVersionString,  // OUT: place reserved for asciiz string
  int pnStringSize);      // IN:  string length

EXTERNC int FFL_FileLength(
  const char* pcFileName);

//== MemBuf handling ===============================================
// callback procedure
// During reading callback is invoked if read requests have exhausted buffer ('hunger').
// The procedure should copy some more data (not more than 'nbytes') to 'buf' pointer,
// and return number of bytes copied there.
// During writing callback is invoked if write requests have overflown buffer ('flood').
// The procedure should read some data (not more than 'nbytes') from 'buf' pointer,
// and return number of bytes it have read from there.
typedef int (*pfnMemIOCallback)(
  void* buf,    // where to read/write data
  int nbytes);  // number of bytes to read/write

DECL_HANDLE(HMemIOBuf);

// create Mem I/O structure and buffer inside it.
// used for writing only.
EXTERNC MIA_RESULT_CODE FFL_BUF_Create(
  HMemIOBuf* psmb,
  uint32 size,
  pfnMemIOCallback* pcb);

// create Mem I/O structure based on given buffer
EXTERNC MIA_RESULT_CODE FFL_BUF_Open(
  HMemIOBuf* psmb,
  void* buffer,
  uint32 size,
  pfnMemIOCallback* pcb,
  int isForWriting);

// freE Mem I/O resources
// Flushes before exit.
EXTERNC MIA_RESULT_CODE FFL_BUF_Close(
  HMemIOBuf psmb);

// flush data gained in Mem I/O structure to external storage
// used for writing only.
EXTERNC MIA_RESULT_CODE FFL_BUF_Flush(
  HMemIOBuf psmb);

// read bytes from Mem I/O
// no argument checks are done to increse performance
EXTERNC MIA_RESULT_CODE FFL_BUF_Read(
  HMemIOBuf psmb,
  void* dst,
  uint32 nb);

// write bytes from Mem I/O
// no argument checks are done to increse performance
EXTERNC MIA_RESULT_CODE FFL_BUF_Write(
  HMemIOBuf psmb,
  const void* dst,
  uint32 nb);

// get size of data currently in Mem I/O structure
EXTERNC int FFL_BUF_GetCurrentSize(
  HMemIOBuf psmb);

// get total number of data passed through Mem I/O (including currently hold data)
EXTERNC int FFL_BUF_GetTotalSize(
  HMemIOBuf psmb);

// determine if there is flood or hunger error
EXTERNC int FFL_BUF_IsFloodOrHunger(
  HMemIOBuf hsmb);

EXTERNC int FFL_BUF_MRead(
  HMemIOBuf hsmb,
  void* dst,
  uint32 nb);

//== BMP handling ==============================================================
#define SINGLEBITPP   0
#define SINGLEBYTEPP  1
#define GRAYSCALE     2
#define RGB24         3
#define BMPTYPEMASK   0xffff
#define MIRROR_X      0x10000
#define MIRROR_Y      0x20000
#define USUAL_BMP     (GRAYSCALE|MIRROR_Y)

EXTERNC MIA_RESULT_CODE FFL_BMP_Save(
  const void* buf,      // pointer to data
  int width,  // size of bitmap in pixels
  int height, //
  int stride, // stride in bytes
  const char* filename, // file name
  int type);  // type of a bitmap / mirroring

// load a bitmap from file
EXTERNC MIA_RESULT_CODE FFL_BMP_Load(
  void* im,             // destination image
  int stride, // stride in bytes
  int  size,  // size of buffer
  const char* filename, // file name
  int type);    // type of a bitmap / mirroring

EXTERNC MIA_RESULT_CODE FFL_BMP_GetType(
  const char* filename,
  int* sizeX,
  int* sizeY,
  int* type);

// read grayscale BMP
EXTERNC MIA_RESULT_CODE FFL_BMP_LoadSimple(
  uint8** pim,        // OUT: allocated buffer with image
  int* pxs,           // OUT: image size
  int* pys,
  const char* name);  // IN:  filename

//== PGM/PPM handling ==============================================================
// get image dimensions from PGM/PPM file
EXTERNC MIA_RESULT_CODE FFL_PGM_GetType(
  uint32 *xs,          // OUT: image size in pixels
  uint32 *ys,          //
  uint32 *nc,          // OUT: number of colors (1 for PGM, 3 for PPM)
  const char* filename);      // IN: name of the file

// save image to PGM/PPM file
EXTERNC MIA_RESULT_CODE FFL_PGM_Save(
  const unsigned char* src,   // IN: source image
  uint32 xs,           // IN: image size in pixels
  uint32 ys,           //
  uint32 nc,           // IN: number of components
  uint32 str,          // IN: image stride in bytes
  const char* filename);      // IN: name of the file

// load image from PGM/PPM file
EXTERNC MIA_RESULT_CODE FFL_PGM_Load(
  unsigned char* dst,     // OUT: destination image
  uint32 str,      // IN:  image buffer stride
  const char* filename);  // IN:  name of the file

//== adaptive Huffman codec =============================================
// encode by self-made adaptive Huffman
EXTERNC int FFL_ADH_Encode(
  unsigned char* obuf,
  int olen,
  const unsigned char* ibuf,
  int ilen);

// decode from self-made adaptive Huffman
EXTERNC int FFL_ADH_Decode(
  unsigned char* obuf,
  int olen,
  const unsigned char* ibuf,
  int ilen);

//== Lempel-Ziv codec ===================================================
// encode by self-made LZW
EXTERNC int FFL_LZW_Encode(
  unsigned char* obuf,
  const unsigned char* ibuf,
  int a_size);

// decode from self-made LZW
EXTERNC int FFL_LZW_Decode(
  unsigned char* obuf,
  const unsigned char* ibuf,
  int* coded);

//== base64 codec =======================================================
// code to base64
#define  FFL_B64_SizeOfEncoded(n) ((((n)/3)*4)+(((n)%3)?(((n)%3)+1):0)+1)

EXTERNC MIA_RESULT_CODE FFL_B64_Encode(
  char* dst,
  int32* dstlen,
  const void* src,
  int32 buflen);

// decode from base64
#define FFL_B64_SizeOfDecoded(n) ((n%4!=2)?(((n-1)/4)*3+((((n-1)%4)+1)>>1)):(-2))

EXTERNC MIA_RESULT_CODE FFL_B64_Decode(
  void* dst,
  int32* dstlen,
  const char* src,
  int32 buflen);

//== zlib codec =======================================================
// compress memory block by a general-pupose zipper
EXTERNC MIA_RESULT_CODE FFL_ZLIB_Encode(
  void* dst,                // IN:     destination buffer
  int32* dst_len,   // IN/OUT: reserved/filled space
  const void* src,          // IN:     source buffer
  int32 src_len);   // IN:     source buffer size

// extract memory block that was compressed general-pupose zipper
EXTERNC MIA_RESULT_CODE FFL_ZLIB_Decode(
  void* dst,                // IN:     destination buffer
  int32* dst_len,   // IN/OUT: reserved/filled space
  const void* src,          // IN:     source buffer
  int32 src_len);   // IN:     source buffer size

//== png codec and file IO ============================================
// compress image by PNG
EXTERNC MIA_RESULT_CODE FFL_PNG_Encode(
  void* dst,                // destination buffer
  int32* dstsize,   // IN: dst buf len, OUT: number of bytes pushed to dst
  const void* src,          // source buffer
  int32 xs,         // image size
  int32 ys,         //
  int32 nc);        // number of components

// expand 8-bit image previously compressed by PNG
EXTERNC MIA_RESULT_CODE FFL_PNG_Decode(
  void* dst,              // IN:  destination buffer
  int32* xs,      // OUT: size of image
  int32* ys,      // OUT:
  int32* nc,      // IN:  dst buf len, OUT: number of color components
  const void* src,        // IN:  source buffer
  int32 src_len); // IN:  source buffer length

EXTERNC MIA_RESULT_CODE FFL_PNG_Load(
  unsigned char* Dest,    // IN:  buffer for image
  int*  pxs,              // OUT: horizontal size
  int*  pys,              // OUT: vertical size
  int*  pnc,              // OUT: number of color components
  int*  pnd,              // IN/OUT: alloced bytes/bit color depth
  const char* FileName);  // IN:  filename

//== JPG codec and file IO =================================================
// compress image by JPG
EXTERNC MIA_RESULT_CODE FFL_JPG_Encode(
  void* outbuffer,              // destination buffer
  int32* dstsize,       // IN: dst buf len, OUT: number of bytes pushed to dst
  const void* image_buffer,     // source buffer
  int32 image_width,    // image size
  int32 image_height,   //
  int32 nc,             // IN: number of color components
  int quality);                 // image compression quality

// expand image previously compressed by JPG
EXTERNC MIA_RESULT_CODE FFL_JPG_Decode(
  void* outbuf,                 // destination buffer
  int32* image_width,   // image size
  int32* image_height,  //
  int32* nc,            // IN:  dst buf len, OUT: number of color components
  const void* inbuf,            // source buffer
  int32 inbuflen_);     // allocated input buffer

EXTERNC MIA_RESULT_CODE FFL_JPG_Load(
  unsigned char* Dest,    // IN:  buffer for image
  int*  pxs,              // OUT: horizontal size
  int*  pys,              // OUT: vertical size
  int*  pnc,              // OUT: number of color components
  int*  pnd,              // IN/OUT: alloced bytes/bit color depth
  const char* FileName);  // IN:  filename

//== JP2K codec and file IO =====================================================
// compress image by JPEG-2000
EXTERNC MIA_RESULT_CODE FFL_JP2_Encode(
  void* dst,              // OUT: destination buffer
  int32* dstsize, // IN: dst buf len, OUT: number of bytes pushed to dst
  const void* src,        // IN: source buffer
  int32 xs,       // IN: image size
  int32 ys,       //
  int32 nc,       // IN: number of color components
  int quality);           // IN: image compression quality, in 1/1000 of initial size

// expand image previously compressed by JPEG-2000
EXTERNC MIA_RESULT_CODE FFL_JP2_Decode(
  void* dst,                // OUT: destination image buffer
  int32* pxs,       // OUT: image size
  int32* pys,       //
  int32* pnc,       // IN/OUT: allocated dst bytes/number of colors
  const void* src,          // IN:  source block
  int32 srcsize);   // IN:  source block size

EXTERNC MIA_RESULT_CODE FFL_JP2_Load(
  unsigned char* Dest,    // IN:  buffer for image
  int*  pxs,              // OUT: horizontal size
  int*  pys,              // OUT: vertical size
  int*  pnc,              // OUT: number of color components
  int*  pnd,              // IN/OUT: alloced bytes/bit color depth
  const char* FileName);  // IN:  filename

//== lossless JPEG codec and file IO ========================================
// compress by LOCO
EXTERNC MIA_RESULT_CODE FFL_JLS_Encode(
  void* dst,                // OUT: destination buffer
  int32* dst_size,  // IN/OUT: allocated/used bytes
  const void* src,          // IN:  source image
  int32 xs,         // IN:  image size
  int32 ys,         //
  int32 nc,         // IN: number of color components
  int32 near_);     // IN: deviation (0 for lossless)

// decompress LOCO block
EXTERNC MIA_RESULT_CODE FFL_JLS_Decode(
  void* dst,                // OUT: destination image buffer
  int32* xs,        // OUT: image size
  int32* ys,        //
  int32* nc,        // IN/OUT: allocated dst bytes/number of colors
  const void* src,          // IN:  source block
  int32 srcsize);   // IN:  source block size

EXTERNC MIA_RESULT_CODE FFL_JLS_Load(
  unsigned char* Dest,    // IN:  buffer for image
  int*  pxs,              // OUT: horizontal size
  int*  pys,              // OUT: vertical size
  int*  pnc,              // OUT: number of color components
  int*  pnd,              // IN/OUT: alloced bytes/bit color depth
  const char* FileName);  // IN:  filename

//== SFALIC codec =====================================================
EXTERNC int FFL_SFL_Encode(
  unsigned char* outbuf,
  int* plen,
  const unsigned char* inbuf,
  int xs,
  int ys);

EXTERNC int FFL_SFL_Decode(
  unsigned char* outbuf,
  int* pxs,
  int* pys,
  const unsigned char* inbuf,
  int len);

//== DTC file handling =========================================================
// file types
#define FFL_DTC_FILETYPE_NONE 0  // no file opened
#define FFL_DTC_FILETYPE_DAT  1       // old-style DAT file (raw data)
#define FFL_DTC_FILETYPE_DTC1 2      // DTC1 file format (raw data+header in the end)
#define FFL_DTC_FILETYPE_DTC2 3      // DTC2 file format header+individually compresssed data frames
#define FFL_DTC_FILETYPE_DTC3 4      // DTC3 2 bytes reserved for frame count instead of 1 in DTC2
#define FFL_DTC_FILETYPE_DTC4 5      // DTC4 additional blocks, 4 bytes for frame count
#define FFL_DTC_FILETYPE_DTC_DEFAULT 5  // default

// compression types
#define FFL_DTC_COMPRESSIONTYPE_NONE 0      // no compression (raw data)
#define FFL_DTC_COMPRESSIONTYPE_ZLIB 1      // zlib'ed
#define FFL_DTC_COMPRESSIONTYPE_PNG  2      // png'ed
//#define FFL_DTC_COMPRESSIONTYPE_TIF  3      // tiff'ed - excluded
#define FFL_DTC_COMPRESSIONTYPE_JPG  4      // jpg'ed
#define FFL_DTC_COMPRESSIONTYPE_JLS  5      // JLS (LOCO)
#define FFL_DTC_COMPRESSIONTYPE_JP2  6      // JPEG-2000 (JPC modification)
#define FFL_DTC_COMPRESSIONTYPE_SFL  7      // SFALIC / (c) R.Starosolski

//++ block domains and data types
#define FFL_DTC_BLOCK_INVALID  0xffffffff
// per-X fields (block domains)
#define FFL_DTC_DOMAIN_IMAGE   0x00000000
#define FFL_DTC_DOMAIN_FRAME   0x10000000
#define FFL_DTC_DOMAIN_CHANNEL 0x20000000
#define FFL_DTC_DOMAIN_FILE    0x40000000
// type of information
#define FFL_DTC_BLOCK_NONE                    0 // dummy
#define FFL_DTC_BLOCK_GENERIC                 1 // just a memory block
#define FFL_DTC_BLOCK_IMAGE                   2 // image
#define FFL_DTC_BLOCK_SOURCEFILENAME          3 // name of source file
#define FFL_DTC_BLOCK_EYECOLOR          0x10001 // color of eye, text
#define FFL_DTC_BLOCK_RACE              0x10002 // race, text
#define FFL_DTC_BLOCK_AGE               0x10003 // age in years, text
#define FFL_DTC_BLOCK_SEX               0x10004 // sex, text
#define FFL_DTC_BLOCK_TIMEMARK          0x20001 // time mark in msec (DM-IriMagic)

// DTC file information structure
typedef struct
{
  int filetype;       // type of file
  // dimension block
  uint16 xs;          // image width in pixels
  uint16 ys;          // image height in pixels
  uint16 nFr;         // number of frames in sequence, char in DTC2, short in DTC3/4
  uint8 nCh;          // number of channels
  uint8 bpP;          // bits per pixel
  uint8 comprtype;    // compression method
  int compression_losses;     // compression quality (JPG, JLS). 0 for lossless
                              // JLS - maximum brightness deviation
                              // JPG - percentage of loss
                              // JP2 - 0.001 parts of original size, 1000 for lossless
  // three 'magic' blocks
  uint32 nCommonLen;   // common parameters length
  char* pCommonParams;        // common parameters (lighter positions, etc.)
  uint32 nFaceLen;     // face parameters length
  char* pFaceParams;          // face parameters
  uint32 nIrisLen;     // iris parameters length
  char* pIrisParams;          // iris parameters
} SDTCFileInfo;

#define FFL_DTC_ClearFileInfo(psdi)  memset(psdi,0,sizeof(SDTCFileInfo))
#define FFL_DTC_FreeFileInfo(psdi)   { if ((psdi)->pCommonParams) free((psdi)->pCommonParams); memset(psdi,0,sizeof(SDTCFileInfo)); };

DECL_HANDLE(HDtcFile);

// create new DTC file (for writing)
EXTERNC MIA_RESULT_CODE FFL_DTC_CreateFile(
  HDtcFile* pDTCHandler,        // OUT: handler to be passed to all subsequent functions
  const char* filename,         // IN:  file name
  const SDTCFileInfo* sdfi);    // IN:  file parameters
// open existing DTC file (for reading only)
EXTERNC MIA_RESULT_CODE FFL_DTC_OpenFile(
  HDtcFile* pDTCHandler,        // OUT: handler to be passed to all subsequent functions
  const char* filename,         // IN:  file name
  SDTCFileInfo* sdfi);          // OUT: (optional) file parameters
// close DTC file
EXTERNC MIA_RESULT_CODE FFL_DTC_CloseFile(
  HDtcFile DTCHandler);         // IN:  handler
// write data to DTC file
EXTERNC MIA_RESULT_CODE FFL_DTC_PutData(
  HDtcFile DTCHandler,   // IN:  handler
  const void* data,      // IN:  data block
  uint32 nCh,            // IN:  channel to put this image to
  uint32 nFr);           // IN:  frame to put this image to
EXTERNC MIA_RESULT_CODE FFL_DTC_PutDataExt(
  HDtcFile DTCHandler,   // IN:  handler
  const void* data,      // IN:  data block
  uint32 len,            // IN:  length of data block (ignored for image type)
  uint32 nCh,            // IN:  channel to put this image to
  uint32 nFr,            // IN:  frame to put this image to
  uint32 nType);         // IN:  type of data
// write data in a form of prepared (compressed) chunk to DTC file directly
EXTERNC MIA_RESULT_CODE FFL_DTC_PutDataDirect(
  HDtcFile DTCHandler,   // IN:  handler
  const void* data,      // IN:  data block
  uint32 nCh,            // IN:  channel to put this image to
  uint32 nFr,            // IN:  frame to put this image to
  uint32 len);           // IN:  chunk length
// read data in a form of prepared (compressed) chunk from DTC file directly
EXTERNC MIA_RESULT_CODE FFL_DTC_GetDataDirect(
  HDtcFile DTCHandler,   // IN:  handler
  void* data,            // OUT: data block
  uint32 nCh,            // IN:  channel to put this image to
  uint32 nFr,            // IN:  frame to put this image to
  uint32* len);          // IN:  reserved chunk length, OUT: really obtained bytes
// read data from DTC file
EXTERNC MIA_RESULT_CODE FFL_DTC_GetData(
  HDtcFile DTCHandler,   // IN:  handler
  void* data,            // IN:  data block
  uint32 nCh,            // IN:  channel to put this image to
  uint32 nFr);           // IN:  frame to put this image to
// read data from DTC file
EXTERNC MIA_RESULT_CODE FFL_DTC_GetDataExt(
  HDtcFile DTCHandler,   // IN:  handler
  void* data,            // OUT: data block
  uint32* plen,          // OUT: length of data block (ignored for image type)
  uint32 nCh,            // IN:  channel to put this image to
  uint32 nFr,            // IN:  frame to put this image to
  uint32 nType);         // IN:  type of data
// read data from DTC file - simple but slow
EXTERNC MIA_RESULT_CODE FFL_DTC_GetDataSimple(
  const char* filename,   // IN:  file name
  void* data,             // IN:  data block
  int32* len,             // IN/OUT:  buffer alloced / BpP
  int32* nCh,             // IN/OUT:  channel/xs
  int32* nFr);            // IN/OUT:  frame/ys
// check presence of frame in DAT/DTC1/DTC2/DTC3/DTC4 file
EXTERNC MIA_RESULT_CODE FFL_DTC_CheckData(
  HDtcFile DTCHandler,    // IN:  handler
  int* bExists,           // OUT: data item exists
  uint32 nCh,             // IN:  channel to put this image to
  uint32 nFr,             // IN:  frame to put this image to
  uint32 type);           // IN:  type of data (for DTC4), ignoreed for others

//== unified (any extension) file reading =========================================
// function for universal reading of BMP/JPG/JLS/JP2/PNG/DTC
// judges file format by extension
EXTERNC MIA_RESULT_CODE FFL_UNI_ReadImage(
  unsigned char* buf,     // IN:  externally alloced buffer for image
  const char* filename,   // IN:  image file name
  int* pxs,               // IN/OUT: channel (for DTC)/x size
  int* pys,               // IN/OUT: frame(for DTC)/y size
  int* pnc);              // IN/OUT: buffer size/BpP

// function for universal reading of BMP/JPG/JLS/JP2/PNG/DTC/DAT
// judges file format by extension, uses augmented names (for DAT/DTC)
MIA_RESULT_CODE FFL_UNI_ReadImage_AugNam(
  unsigned char* buf,     // IN:  externally alloced buffer for image
  const char* filename,   // IN:  image file name (augmented in DTC case)
  int* pxs,               // OUT: x size
  int* pys,               // OUT: y size
  int* pnc);              // IN/OUT: buffer size/BpP

//== DT0 file handling =========================================================
// DT0 file information structure
typedef struct
{
  FILE* fp;                 // file pointer
  unsigned int isWritten;   // is file created for writing
  unsigned int xs;          // image width in pixels
  unsigned int ys;          // image height in pixels
  unsigned int nFr;         // number of frames in sequence
  unsigned int cFr;         //
  unsigned int isLeft;      // =0 for right and =1 for left frame
} SDT0FileInfo;

#define SDT0FileInfo_Clear(psdi)  MIA_memset(psdi,0,sizeof(SDTCFileInfo))

// create new DT0 file
EXTERNC MIA_RESULT_CODE FFL_DT0_CreateFile(
  SDT0FileInfo* sdfi,           // IN/OUT: structure
  const char* filename);        // IN:     file name

// open existing DT0 file
EXTERNC MIA_RESULT_CODE FFL_DT0_OpenFile(
  SDT0FileInfo* sdfi,           // OUT: structure
  const char* filename);        // IN:  file name

// close DT0 file
EXTERNC MIA_RESULT_CODE FFL_DT0_CloseFile(
  SDT0FileInfo* sdfi);          // IN:  handler

// write data to DT0 file
EXTERNC MIA_RESULT_CODE FFL_DT0_PutData(
  SDT0FileInfo* sdfi,           // IN:  handler
  const void* data,             // IN:  data block
  unsigned int nFr);            // IN:  frame to put this image to

// read data from DTC file
EXTERNC MIA_RESULT_CODE FFL_DT0_GetData(
  SDT0FileInfo* sdfi,           // IN:  handler
  void* data,                   // IN:  data block
  unsigned int nFr);            // IN:  frame to put this image to

//== IIIF file handling =================================================

//-- constants -------------
// capture devices
#define IIIF_CAPTURE_DEVICE_UNDEF                      0

// image properties (bitfield)
#define IIIF_HORIZONTAL_ORIENTATION_UNDEF         0x0000
#define IIIF_HORIZONTAL_ORIENTATION_BASE          0x0001
#define IIIF_HORIZONTAL_ORIENTATION_FLIPPED       0x0002
#define IIIF_VERTICAL_ORIENTATION_UNDEF           0x0000
#define IIIF_VERTICAL_ORIENTATION_BASE            0x0004
#define IIIF_VERTICAL_ORIENTATION_FLIPPED         0x0008
#define IIIF_SCAN_TYPE_UNDEF                      0x0000
#define IIIF_SCAN_TYPE_PROGRESSIVE                0x0010
#define IIIF_SCAN_TYPE_INTERLACE_FRAME            0x0020
#define IIIF_SCAN_TYPE_INTERLACE_FIELD            0x0030
#define IIIF_SCAN_TYPE_CORRECTED                  0x0040
#define IIIF_IROCC_UNDEF                          0x0000
#define IIIF_IROCC_ZEROFILL                       0x0080
#define IIIF_IROCC_UNITFILL                       0x0180
#define IIIF_IRBNDY_UNDEF                         0x0000
#define IIIF_IRBNDY_PROCESSED                     0x0200

// image formats
#define IIIF_IMAGEFORMAT_MONO_RAW                      2
#define IIIF_IMAGEFORMAT_RGB_RAW                       4
#define IIIF_IMAGEFORMAT_MONO_JPEG                     6
#define IIIF_IMAGEFORMAT_RGB_JPEG                      8
#define IIIF_IMAGEFORMAT_MONO_JPEG_LS                 10
#define IIIF_IMAGEFORMAT_RGB_JPEG_LS                  12
#define IIIF_IMAGEFORMAT_MONO_JPEG2000                14
#define IIIF_IMAGEFORMAT_RGB_JPEG2000                 16

// image sizes
#define IIIF_WIDTH_UNDEF                               0
#define IIIF_HEIGHT_UNDEF                              0
#define IIIF_INTENSITY_DEPTH_UNDEF                     0

// transformation to polar image algorithm
#define IIIF_TRANS_UNDEF                               0
#define IIIF_TRANS_STD                                 1

// feature identifiers
#define IIIF_EYE_UNDEF                                 0
#define IIIF_EYE_RIGHT                                 1
#define IIIF_EYE_LEFT                                  2

// image quality
#define IIIF_IMAGE_QUAL_UNDEF                       0xfe
#define IIIF_IMAGE_QUAL_LOW                           26
#define IIIF_IMAGE_QUAL_MED                           51
#define IIIF_IMAGE_QUAL_HIGH                          76
#define IIIF_IMAGE_QUAL_BEST                         100

// rotations
#define IIIF_ROT_ANGLE_UNDEF                      0xffff
#define IIIF_ROT_UNCERTAIN_UNDEF                  0xffff

//-- data exchange structures--//
// image info
typedef struct
{
  unsigned char   nQuality;                       // quality of this image
  unsigned short  nRotationAngle;                 // angle of rotation
  unsigned short  nRotationUncertainty;           // uncertainty of rotation
  // inner data
  void*           pvReserved;                     // user should never touch this
} SIIIFImageInfo;

// feature info
typedef struct
{
  unsigned char   nFeatureIdentifier;             // right, left or undef
  unsigned short  nImagesNumber;                  // number of images for this feature
  // per-image data
  SIIIFImageInfo* psImgInfo;                      // info's for each image
} SIIIFeatureInfo;

// record info (for opening)
typedef struct
{
  uint32   nCbeffID;                       // CBEFF ID of vendor
  unsigned short  nCaptureDeviceIdentifier;       // vendor ID of device
  unsigned char   nFeaturesNumber;                // number of features i.e. eyes
  unsigned short  nImageProperties;               // orientation, scan type, etc.
  unsigned short  nIrisDiameter;                  // expected iris diameter in pixels
  unsigned short  nImageFormat;                   // image format
  unsigned short  nRawImageWidth;                 // image size
  unsigned short  nRawImageHeight;                //
  unsigned char   nIntensityDepth;                // bits per color
  unsigned char   nImageTransformation;           // transformation to polar type
  unsigned char   aucDeviceUniqueIdentifier[16];  // DUID
  unsigned char   aucGlobalUniqueIdentifier[16];  // GUID
  // per-feature data
  SIIIFeatureInfo* psFeatInfo;                    // info's for each feature
} SIIIFRecordInfo;

typedef void* HIIIFHandle;

//-- functions ------------------------
// create new IIIF file (for writing)
EXTERNC MIA_RESULT_CODE FFL_IIIF_CreateFile(
  HIIIFHandle* pHIIIF,          // OUT: handler to be passed to all subsequent functions
  const char* pcFilename,       // IN:  file name
  const SIIIFRecordInfo* pSIFR);// IN:  file parameters
// open existing IIIF file (for reading only)
EXTERNC MIA_RESULT_CODE FFL_IIIF_OpenFile(
  HIIIFHandle* pHIIIF,          // OUT: handler to be passed to all subsequent functions
  const char* pcFilename,       // IN:  file name
  SIIIFRecordInfo** ppSIFR);    // OUT: file parameters
// close IIIF file
EXTERNC MIA_RESULT_CODE FFL_IIIF_CloseFile(
  HIIIFHandle HIIIF);           // IN:  handler
// add feature description to iiif file
EXTERNC MIA_RESULT_CODE FFL_IIIF_AddFeature(
  HIIIFHandle HIIIF,                // IN:  handler
  const SIIIFeatureInfo* pFeatInfo, // IN:  feature info
  unsigned int nFeature);           // IN:  number of feature in the list
// write data to IIIF file
EXTERNC MIA_RESULT_CODE FFL_IIIF_AddImage(
  HIIIFHandle HIIIF,            // IN:  handler
  const unsigned char* pucImage,// IN:  image data block
  const SIIIFImageInfo* pSIFI,  // IN:  image info
  unsigned int nCh,             // IN:  channel (i.e. feature) to put this image to
  unsigned int nFr);            // IN:  frame to put this image to
// read image from IIIF file
EXTERNC MIA_RESULT_CODE FFL_IIIF_GetImage(
  HIIIFHandle HIIIF,            // IN:  handler
  unsigned char* pucImage,      // IN:  image data block
  unsigned int nCh,             // IN:  channel to put this image to
  unsigned int nFr);            // IN:  frame to put this image to

// load data from TIF file
EXTERNC MIA_RESULT_CODE FFL_TIF_Load(
  unsigned char** ppIm, // OUT: image data, can be NULL if data not needed
  int* pxs,             // OUT: width, can be NULL if data not needed
  int* pys,             // OUT: height, can be NULL if data not needed
  int* pnc,             // OUT: number of color components, can be NULL if no need
  int* pnd,             // OUT: bits per color component, can be NULL if no need
  const char* filename);// IN:  file name

// release the buffer allocated inside FFL library, set pointer to NULL
EXTERNC MIA_RESULT_CODE FFL_Free_Buffer(
  void** ppvBuf);

// recursive scanning and compression of files fitting to masks
EXTERNC MIA_RESULT_CODE FFL_ZIM_CompressFiles(
  const char* archname,     // name of target archive file
  const char* dir,          // root source directory
  const char** masks,       // list of masks. end of list signed by NULL pointer
  time_t beg_tim,           // beginning time
  time_t end_tim,           // ending time
  uint32 SrcChunkSize,      // preferred size of chunk for compression
  int doRecusrion);         // if directory recursion should be done

// restoring the compressed structure
EXTERNC MIA_RESULT_CODE FFL_ZIM_RestoreFiles(
  const char* dir_,         // name of target root directory
  const char* archname);    // name of target source file

#endif //__ffl_h__
