#include <math.h>
#include <stdio.h>
#include "stddefs.h"
#include "Geometry.h"

/*
typedef struct CErr
{
	double sum, ctr, pup, iri;
} CErr;
*/

typedef struct SError
{
	double total, iris, pupil, pCenter, iCenter;
	int hSum[100] = { 0 };
	int hcP[100] = { 0 };
	int hP[100] = { 0 };
	int hcI[100] = { 0 };
	int hI[100] = { 0 };
};

typedef struct SMarkup
{
	char filename[MAX_FILENAME];
	SCircleData pupil, iris;
} SMarkup;


int getExpertMarkup(FILE* fin, CInfo* pup, CInfo* iri, char* name, char* out);
int calcError(CErr* err, CInfo* myPup, CInfo* myIri, CInfo* pup, CInfo* iri);