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

//EXTERNC int ibd_grapar(CInfo* CI, char* name, const unsigned char* img, int H, int W, int dPhi, SBoundPnt* pnBoundary, CInfo* CP, int flags);

EXTERNC int ibd_IrisSegmentation(CInfo* CI,CInfo* CP,CInfo *rCP, char* name, unsigned char* imgInput,  unsigned char *imgOutput, int H, int W, int angle, FILE* log, int flags);

int FLT_Sobel3x3(int* dest,//destination image
				 unsigned char* imgSobel,//gradient magnitude image
				 const unsigned char *img,//original image
				 char *name,//original image name
				 int H,//height 
				 int W,//width
				 int* mask,	//Sobel mask using for gradient calculation
				 int* imdx,	//gradient values in x direction
				 int* imdy,	//gradient values in y direction
				 int flags);

EXTERNC int FLT_Sobel10x10(int* dest,unsigned char* imgSobel, const unsigned char *img, char *name, int H, int W, int* imdx, int* imdy, int flags);
EXTERNC int ibd_graparBothBorders(CInfo* CI, CInfo* CP, char* name, const unsigned char* img, unsigned char *imgBl, int H, int W, int dPhi, int flags);

EXTERNC int IBD_RefinePupil(int* xpList, int *ypList, int* radmean, CInfo* CP, int MinRad, int MaxRad, char *name, const unsigned char* img, int H, int W, int flags);
EXTERNC int IBD_RefineIris(CInfo *CP, CInfo *CI, const unsigned char *img, char* name, int H, int W, int flags);


