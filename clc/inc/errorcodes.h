/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  14 March 2001                                          */
/*      Modified: many times...                                          */
/*      Last revision: 8 February 2015 - file number increased to 4095   */
/*                                       error number decreased to 255   */
/*      Purpose:  Error report codes used everywhere in MIA soft         */
/*      Authors:                                                         */
/*        Ivan Matveev                                                   */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#ifndef __errorcodes_h__
#define __errorcodes_h__

/*
Format of error code is:
Bits: 33222222222211111111110000000000
Bits: 10987654321098765432109876543210
      --------------------------------
      ffffffffffffllllllllllllnnnnnnnn
where f - file number field (up to 4095 files)
      l - line number field (up to 4095 lines)
      n - error number (up to 256 errors)
*/

#ifndef __FILENUM__
  #define __FILENUM__ 0 // default (__FILENUM__TAG0)
#endif

#define ENHANCED_ERROR_CODE(v) \
  ( (((__FILENUM__)&0xfff)<<20) + ((__LINE__&0xfff)<<8) + ((v)&0xff) )

/* common error codes -------------------------------------------------*/
/* everything OK */
#define ERR_OK 0
/* function not implemented (develop-time usage) */
#define ERR_GEN_NOTIMPL ENHANCED_ERROR_CODE(0x01)
/* internal error (problem in library, not callers fault) */
#define ERR_GEN_INTERNAL ENHANCED_ERROR_CODE(0x02)
/* NULL pointer unexpected */
#define ERR_GEN_NULLPOINTER ENHANCED_ERROR_CODE(0x03)
/* data cannot fit in buffer */
#define ERR_GEN_INSUFFICIENT_BUFFER ENHANCED_ERROR_CODE(0x04)
/* sizes do not match */
#define ERR_GEN_SIZE_NOT_MATCH ENHANCED_ERROR_CODE(0x05)
/* bad data alignment */
#define ERR_GEN_BAD_ALIGNMENT ENHANCED_ERROR_CODE(0x06)
/* bad object ID detected */
#define ERR_GEN_INVALID_OBJECT_ID ENHANCED_ERROR_CODE(0x07)
/* bad parameter passed to function */
#define ERR_GEN_INVALID_PARAMETER ENHANCED_ERROR_CODE(0x08)
/* object/variable has not been initialized */
#define ERR_GEN_NOT_INITIALISED ENHANCED_ERROR_CODE(0x09)
/* object/variable is initialized already */
#define ERR_GEN_INITIALISED_ALREADY ENHANCED_ERROR_CODE(0x0a)
/* nothing to proceed (ROI empty, etc) */
#define ERR_GEN_NO_DATA ENHANCED_ERROR_CODE(0x0b)
/* no memory */
#define ERR_GEN_NOMEMORY ENHANCED_ERROR_CODE(0x0c)
/* bad data type */
#define ERR_GEN_BADDATATYPE ENHANCED_ERROR_CODE(0x0d)
/* empty ROI specified */
#define ERR_GEN_EMPTY_ROI ENHANCED_ERROR_CODE(0x0e)
/* error in data */
#define ERR_GEN_BAD_DATA ENHANCED_ERROR_CODE(0x0f)

/* Database errors ---------------------------------------------------*/
/* cannot create (or open for writing) */
#define ERR_DBE_CANNOT_CREATE_FILE ENHANCED_ERROR_CODE(0xd0)
/* cannot open existing file (for reading) */
#define ERR_DBE_CANNOT_OPEN_FILE ENHANCED_ERROR_CODE(0xd1)
/* cannot write to file */
#define ERR_DBE_CANNOT_WRITE_FILE ENHANCED_ERROR_CODE(0xd2)
/* cannot read from file */
#define ERR_DBE_CANNOT_READ_FILE ENHANCED_ERROR_CODE(0xd3)
/* similar item already exists */
#define ERR_DBE_SIMILAR_ITEM_EXISTS ENHANCED_ERROR_CODE(0xd4)
/* cannot find item */
#define ERR_DBE_CANNOT_FIND_ITEM ENHANCED_ERROR_CODE(0xd5)
/* compression error */
#define ERR_DBE_COMPRESSION_ERROR ENHANCED_ERROR_CODE(0xd6)
/* decompression error */
#define ERR_DBE_DECOMPRESSION_ERROR ENHANCED_ERROR_CODE(0xd7)

#define MIA_RESULT_CODE int

#endif /* __errorcodes_h__ */
