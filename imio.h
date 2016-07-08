#include <windows.h>
#include <conio.h>
#define MAX_FILENAME 256

typedef struct
{
	double sum, ctr, pup, iri;
} CErr;

EXTERNC int readBmp8(FILE* imgFile, unsigned char** img,int* H, int* W);
EXTERNC int getExpertMarkup(FILE* fin, CInfo* pup, CInfo* iri, char* name,char* out);
EXTERNC int calcError(CErr* err, CInfo* myPup, CInfo* myIri, CInfo* pup, CInfo* iri);


EXTERNC int CreateBmp8 (char *fname, int Width,int Height, unsigned char* map, BYTE color);
EXTERNC void DRAW_Line(int* mat,int H, int W, int x1, int y1, int x2, int y2);
EXTERNC int SaveBmp8 (char *fname, char* label, int Width,int Height, unsigned char* map, BYTE color);


// draw the circle of given color in grayscale image
EXTERNC void DRAW_CircleInGray(
  unsigned char* im,      // image
  int W,       // image size
  int H,       // image size
  int xc,       // circle center
  int yc,       // circle center
  int r,        // circle radius
  unsigned char color);    // color

EXTERNC void DRAW_SequenceInGray(
	unsigned char* im,		// image
	int W,		//image size
	int H,		//image size
	int* xP,		//x coordinates for the sequence
	int* yP,		//y coordinates for the sequence
	int N,			//Number of elements in the sequence
	unsigned char color);	//color

EXTERNC void Dilate3x3Cross(unsigned char* dst, const unsigned char* src, int W, int H);