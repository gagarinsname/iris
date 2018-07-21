#include "errorcodes.h"
#include "ipl.h"
#include "unitypes.h"

#define BDL_PGRAD_CALCRADIUS    0x00000001  // perform radius calculation
#define BDL_PGRAD_SAVEACC		0x00000002	// save the accumulator image
#define BDL_PGRAD_SAVECANNY		0x00000004	// save the canny edge image
#define BDL_PGRAD_LINEVOTE		0x00000010	// voting in the accumulator using lines
#define BDL_PUPREF_SAVECROP		0x00000020	//save the cropped border image for the pupil
#define BDL_PUPREF_SAVEPOLAR	0x00000040	//save the polar transformed pupil border region
#define FLT_SOBEL_SAVEIMAGE		0x00000080	//save the result of filter applying as an image
#define FLT_SOBEL_SELECTION		0x00000100	//select pixels with the top 15% gradient magnitudes



typedef struct
{
  int xc;   // X coordinate
  int yc;   // Y coordinate
  int r;    // radius (elliptic big axis)
} CInfo;

typedef struct
{
  int16 x;
  int16 y;
  int val;
} SBoundPnt;