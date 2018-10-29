#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <windows.h>
#include <conio.h>
#include <time.h>
#include "header.h"
#include "errorcodes.h"
#include "ipl.h"
#include "unitypes.h"
#include "imio.h"
#include "imPlot.h"
#include "PairedGradients.h"
#include "bdl.h"
#include "Experiment.h"

#define INPUT_FOLDER	"D:/data/"
#define OUTPUT_FOLDER	"D:/data/res/0619/"

#define DB_CASIA2		"params_casia2.txt"
#define DB_IITDELHI		"params_iitdelhi.txt"
#define DB_MIXED		"params_res.txt"

//#include "bpl.h"
//#include "bdl.h"
//#include "dbg.h"

//unsigned char palitra[256][4];

int getMarkupData(SMarkup* cMarkupData, FILE* markupFile)
{
	int i, res;
	char* tmp = (char*)malloc(10 * MAX_FILENAME);
	if (!fgets(tmp, MAX_FILENAME, markupFile))
	{
		fprintf(stderr, "[ERROR] getMarkupData() cannot read another line from markup file.");
		return -1;
	}

	fscanf(markupFile, "%s", tmp);
	strcpy(cMarkupData->filename, tmp);
	
	// skipping unnecessary input
	for (i = 0; i < 11; i++)
		fscanf(markupFile, "%s", tmp);

	fscanf(markupFile, "%s", tmp);
	cMarkupData->pupil.xc = atoi(tmp);
	fscanf(markupFile, "%s", tmp);
	cMarkupData->pupil.yc = atoi(tmp);
	fscanf(markupFile, "%s", tmp);
	cMarkupData->pupil.r = atoi(tmp);
	
	fscanf(markupFile, "%s", tmp);

	fscanf(markupFile, "%s", tmp);
	cMarkupData->iris.xc = atoi(tmp);
	fscanf(markupFile, "%s", tmp);
	cMarkupData->iris.yc = atoi(tmp);
	fscanf(markupFile, "%s", tmp);
	cMarkupData->iris.r = atoi(tmp);
	
	//	printf("Expert:\n(%d ; %d) r = %d\n(%d ; %d) r = %d\n", pup->xc, pup->yc, pup->r, iri->xc, iri->yc, iri->r);
	free(tmp);
	return 0;
}

int calculateError(SError *error, SSegmentationResult *result, SMarkup *markup)
{
	if (result->IrisData->sCircle->xc != -1 && result->IrisData->sCircle->yc != -1)
	{
		error->iCenter = sqrt((double)((result->IrisData->sCircle->xc - markup->iris.xc) * (result->IrisData->sCircle->xc - markup->iris.xc) +
			(result->IrisData->sCircle->yc - markup->iris.yc) * (result->IrisData->sCircle->yc - markup->iris.yc)));
		error->iCenter /= markup->iris.r;
		error->iris = error->iCenter + (double)(result->IrisData->sCircle->r - markup->iris.r) / markup->iris.r;
	}
	else
	{
		error->iCenter = -1.0;
		error->iris = -1.0;
		error->total = -1.0;
	}
	if (result->PupilData->sCircle->xc != -1 && result->PupilData->sCircle->yc != -1)
	{
		error->pCenter = sqrt((double)((result->PupilData->sCircle->xc - markup->pupil.xc) * (result->PupilData->sCircle->xc - markup->pupil.xc) +
			(result->PupilData->sCircle->yc - markup->pupil.yc) * (result->PupilData->sCircle->yc - markup->pupil.yc)));
		error->pCenter /= markup->iris.r;
		error->pupil = error->pCenter + (double)(result->PupilData->sCircle->r - markup->pupil.r) / markup->iris.r;
	}
	else
	{
		error->pCenter = -1.0;
		error->pupil = -1.0;
		error->total = -1.0;
	}

	if (result->PupilData->sRCircle->xc != -1 && result->PupilData->sRCircle->yc != -1)
	{
		error->pCenter = sqrt((double)((result->PupilData->sRCircle->xc - markup->pupil.xc) * (result->PupilData->sRCircle->xc - markup->pupil.xc) +
			(result->PupilData->sRCircle->yc - markup->pupil.yc) * (result->PupilData->sRCircle->yc - markup->pupil.yc)));
		error->pCenter /= markup->iris.r;
		error->pupil = error->pCenter + (double)(result->PupilData->sRCircle->r - markup->pupil.r) / markup->iris.r;
	}

	if (error->pupil != -1 && error->iris != -1)
		error->total = error->pupil + error->iris;

	return 0;
}


int processMarkup(
	int imgBeg,		// a line to start with
	int imgEnd,		// a line to end
	FILE* fout		// to write log file
)
{
	int res = 0;
	char *markupFileName = (char*)malloc(MAX_FILENAME * sizeof(char));
	char *imagePath;
	char *imageOutputPath;
	FILE *imageFile, *markupFile;
	uint8 *imageData = NULL;

	strcpy(markupFileName, INPUT_FOLDER);
	strcat(markupFileName, DB_MIXED);

	int cErrorHist[100] = { 0 };

	if ((markupFile = fopen(markupFileName, "r"))== NULL)
	{
		fprintf(stderr, "[ERROR] processMarkup() wrong filename for markup file.\n");
		return 1;
	}
	fgets(markupFileName, MAX_FILENAME, markupFile);
	imagePath = (char*)malloc(MAX_FILENAME * sizeof(char));
	imageOutputPath = (char*)malloc(MAX_FILENAME * sizeof(char));

	for (int i = 0; i < imgEnd; ++i)
	{
		if (imgBeg > imgEnd)
			break;

		SMarkup cMarkupData;
		SSegmentationResult Result;
		SError error;

		memset(imagePath, 0, MAX_FILENAME * sizeof(char));
		strcpy(imagePath, INPUT_FOLDER);
		memset(imageOutputPath, 0, MAX_FILENAME * sizeof(char));
		strcpy(imageOutputPath, OUTPUT_FOLDER);

		if ((res = getMarkupData(&cMarkupData, markupFile)) != 0)
		{
			fprintf(stderr, "[ERROR]: getMarkupData markup extraction fail.\n");
		}

		if (i < imgBeg)
			continue;
		
		strcat(imagePath, cMarkupData.filename);
		strcat(imageOutputPath, cMarkupData.filename);

		if ((imageFile = fopen(imagePath, "rb")) == NULL)
		{
			fprintf(stderr, "[ERROR] Cannot open the image file.\n");
			res = 3;
		}

		printf("Markup entry %d: %s\n", i, cMarkupData.filename);
		//Read the image data
		int H = -1, W = -1;
		if ((res = readBmp8(imageFile, &imageData, &H, &W)) != 0)
		{
			printf("Error: Can't read the image data.\n");
			res = 4;
		}
		fclose(imageFile);
		if (res != 0)
			break;

		//Iris segmentation function
		if ((res = IS_Initialize(&Result, imageOutputPath)) != 0) {
			printf("\nError: iris segmentation init failed.");
			res = 5;
			break;
		}

		int angle = 120;
		int flags = BDL_PGRAD_CALCRADIUS;
		if ((res = IrisSegmentation(&Result, imageData, H, W, angle, flags)) != 0)
		{
			printf("\nError: Iris segmentation fail.\n");
			res = 6;
			break;
		}

		res = calculateError(&error, &Result, &cMarkupData);
		printf("\nSegmentation error: pCenter = %f, iCenter = %f.\n\n", error.pCenter, error.iCenter);
		
		if ((int)(error.pCenter * 100) < 100)
			++cErrorHist[((int)(error.pCenter * 100))];

		if (error.iCenter > 0.)
		{
			DRAW_CircleInGray(imageData, W, H, Result.IrisData->sCircle, 255);
			DRAW_2DSequenceInGray(imageData, W, H, Result.PupilData->sContour, (const int)(CONTOUR_QUANTIZATION), 255);
			//DRAW_CircleInGray(imageData, W, H, Result.PupilData->sRCircle, 255);
			//DRAW_CircleInGray(imageData, W, H, Result.PupilData->sCircle, 200);
			SaveBmp8(imageOutputPath, "_CIRCLE", W, H, imageData, 4);
		}

		if ((res = IS_Deinitialize(&Result)) != 0)
		{
			printf("\nError: Iris segmentation deinit failed.");
			res = 7;
			break;
		}
	}

	fprintf(fout, "Center error detection histogram.");
	for (int h = 0; h < 100; ++h)
		fprintf(fout, "%d ", cErrorHist[h]);
	fprintf(fout, "\n");

	if (imageData != NULL)
		free(imageData);
	free(imagePath);
	free(imageOutputPath);
	free(markupFileName);
-
	fclose(markupFile);
	return res;
}

int main(int argc, char** argv)
{
	// load markup file
	int res;
	char logFile[MAX_FILENAME] = { 0 };
	strcpy(logFile, OUTPUT_FOLDER);
	FILE* fout = fopen(strcat(logFile, "log.txt"), "wb");
	if ((res = processMarkup(2042, 2050, fout)) != 0)
	{
		fprintf(stderr, "[ERROR] processMarkup returned %d.", res);
	}
	fclose(fout);
	/*
	load markup file
	process markup file
	collect statistics
	*/
	return 0;
}

