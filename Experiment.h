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
};

typedef struct SMarkup
{
	char filename[MAX_FILENAME];
	SCircleData pupil, iris;
} SMarkup;


EXTERNC int getExpertMarkup(FILE* fin, CInfo* pup, CInfo* iri, char* name, char* out);
EXTERNC int calcError(CErr* err, CInfo* myPup, CInfo* myIri, CInfo* pup, CInfo* iri);


