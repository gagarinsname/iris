#include <windows.h>
#include <conio.h>
#include "unitypes.h"
#include "Geometry.h"
#include "GlobalParams.h"

#define MAX_FILENAME 256

typedef struct CErr
{
	double sum, ctr, pup, iri;
} CErr;


EXTERNC int readBmp8(FILE* imgFile, unsigned char** img,int* H, int* W);
EXTERNC int getExpertMarkup(FILE* fin, CInfo* pup, CInfo* iri, char* name,char* out);
EXTERNC int calcError(CErr* err, CInfo* myPup, CInfo* myIri, CInfo* pup, CInfo* iri);
EXTERNC int CreateBmp8 (char *fname, int Width,int Height, unsigned char* map, BYTE color);
EXTERNC int SaveBmp8 (char *fname, char* label, int Width,int Height, unsigned char* map, BYTE color);