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
#include "IrisSegmentation.h"
#include "bdl.h"

#define INPUT_FOLDER "D:/data/casia2/"
#define OUTPUT_FOLDER "D:/data/res/0224/"

//#include "bpl.h"
//#include "bdl.h"
//#include "dbg.h"

unsigned char palitra[256][4];

int test_segm(int argc, char** argv)
{
	clock_t time;
	unsigned char *imgInput, *imgOutput;
	char *imgName, *imgOut, *name, *cname;
	char tmp[MAX_FILENAME];
	int i, k, res, H, W, angle, flags;

	FILE *imgFile;

	name = (char*)malloc(1024 * sizeof(char));
	strcpy(name, OUTPUT_FOLDER);


	imgName = (char*)malloc(256 * sizeof(char));
	strcpy(imgName, INPUT_FOLDER);
	imgOut = (char*)malloc(256 * sizeof(char));
	strcpy(imgOut, OUTPUT_FOLDER);
	
	strcat(imgName, "1000(R-1-0)041108-222400_e_00_0000.bmp");
	strcat(imgOut, "1000(R-1-0)041108-222400_e_00_0000.bmp");

	imgFile = fopen(imgName, "rb");
	if (imgFile == NULL)
	{
		printf("Error: Cannot open the image file.\n");
		return -1;
	}

	//Read the image data
	if ((res = readBmp8(imgFile, &imgInput, &H, &W)) != 0)
	{
		printf("Error: Can't read the image data.\n");
		return -2;
	}
	imgOutput = (uint8*)malloc(H*W * sizeof(uint8));
	memcpy(imgOutput, imgInput, H*W * sizeof(uint8));
	
	//Reduce the file extension
	memcpy(name, imgOut, strlen(imgOut) - 4);
	name[strlen(imgOut) - 4] = 0;
	
	time = clock();

	//Iris segmentation function
	SSegmentationResult Result;// = (SSegmentationResult*)malloc(sizeof(SSegmentationResult));
	if ((res = IS_Initialize(&Result)) != 0) {
		printf("\nError: iris segmentation init failed.");
		return -1;
	}

	angle = 120;
	flags = BDL_PGRAD_CALCRADIUS;
	if ((res = IrisSegmentation(&Result, imgInput, H, W, angle, flags)) != 0)
	{
		printf("\nError: Iris segmentation fail.\n");
		return -1;
	}
	time = clock() - time;
	//getStatistics(SSegmentationResult* Result);

	DRAW_CircleInGray(imgOutput, W, H, Result.IrisData->sCircle, 255);
	//free(Result);
	
	if ((res = IS_Deinitialize(&Result)) != 0)
	{
		printf("\nError: Iris segmentation deinit failed.");
		return -1;
	}

	printf("name: %s", name);
	//SaveBmp8(name, "_pupil_refined.bmp", W, H, imgOutput, 4);

	printf("\ntime: %f seconds.\n\n", (double)time / CLOCKS_PER_SEC);

	free(imgOut);
	free(imgOutput);
	free(imgInput);
	free(imgName);
	fclose(imgFile);
	free(name);
	//fclose(out);
	return 0;
}

