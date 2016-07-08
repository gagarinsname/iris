#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "header.h"
#include "imio.h"
#include "ipl.h"

static int sort_int(const void* e1, const void* e2)
{
    return ((int*)e1)-((int*)e2);
}

static int sort_long(const void* e1, const void* e2)
{
    if (*(long int*)e1 == *(long int*)e2)
		return 0;
	if (*(long int*)e1 > *(long int*)e2)
		return 1;
	else
		return -1;

}

int FLT_Sobel10x10(long* dest, unsigned char* imgSobel, const unsigned char *img, char *name, int H, int W, long* imdx, long* imdy, int flags)
{
	FILE* log = fopen("../data/res/0629/sobel_log.txt", "wb");
	int pad,size,x,y, cPt;
	long int lowMag, highMag, max;
	double val;
	int* imgPad;
	long int *imgHist;
	memset(dest, 0, H*W*sizeof(unsigned char));
	memset(imdx, 0, H*W*sizeof(long));
	memset(imdy, 0, H*W*sizeof(long));
	memset(imgSobel,(unsigned char)0, H*W*sizeof(unsigned  char));
	//Padding

	pad = 11;
	imgPad = (unsigned int*)malloc((H+2*pad)*(W+2*pad)*sizeof(int));
	for (y = 0; y < H; y++)
	{
		for (x = 0; x < W; x++)
		{
			imgPad[(y+pad) * W + (x+pad)] = (int)img[y * W + x];
		}
	}
	for (y = 0; y < pad; y++)
	{
		for (x = pad; x < W + pad; x++)
		{
			imgPad[y * W + x] = (int)img[(2 * pad - y) * W + x];
			//
			imgPad[(H+y) * W + x] = (int)img[(H-y) * W + x];
			//
		}
	}
	for (x = 0; x < pad; x++)
	{
		for (y = pad; y < H+pad; y++)
		{
			imgPad[y * W + x] = (int)img[y * W + (2 * pad-x)];
			imgPad[y * W + W + x] = (int)img[y*W + W - x];
		}
	}
	for (y = 0; y < pad; y++)
	{
		for (x = 0; x < pad; x++)
		{
			imgPad[y * W + x] = (imgPad[y * W + (2*pad-x)] + imgPad[(2*pad - y) * W + x]) / 2;
			imgPad[(H+y) * W + x] = (imgPad[(H+y) * W + (2*pad-x)] + imgPad[(H-y) * W + x]) / 2;
			imgPad[y * W + W + x] = (imgPad[y * W + (W-x)] + imgPad[(pad - y) * W + x+W]) / 2;
			imgPad[(H+y) * W + W + x] = (imgPad[(H+y) * W + (W-x)] + imgPad[(H-y) * W + x+W]) / 2;
		}
	}
	//Gradient X and Y calculation

	/*
	Using big sobel mask:
	1  1  2  2  4  5  7  8  10   11  11  11  10  8  7  5  4  2  2  1  1
	0  0  0  0  0  0  0  0   0    0   0   0   0  0  0  0  0  0  0  0  0 
	0 0 ............................................................. 0
	...................................................................
	0  0  0  0  0  0  0  0   0    0   0   0   0  0  0  0  0  0  0  0  0 
	-1 -1 -2 -2 -4 -5 -7 -8 -10 -11 -11 -11 -10 -8 -7 -5 -4 -2 -2 -1 -1
	*/
	/*
	for (y = 0; y < 2*pad; y++)
	{
		for (x = 0; x < 2*pad; x++)
			//fprintf(log, "%d ", imgPad[y*W+x]);
		//fprintf(log, "\n");
	}*/
	size = 10;
	max = 0;
	
	for (y = 0; y < H; y++)
	{
		for (x = 0; x < W; x++)
		{
			imdy[y*W+x] = (long int)(1*imgPad[(y+pad-size)*W+(x+pad-10)] + 1*imgPad[(y+pad-size)*W+(x+pad-9)] + 2*imgPad[(y+pad-size)*W+(x+pad-8)] + 2*imgPad[(y+pad-size)*W+(x+pad-7)] +
				4*imgPad[(y+pad-size)*W+(x+pad-6)] + 5 * imgPad[(y+pad-size)*W+(x+pad-5)] + 7 * imgPad[(y+pad-size)*W+(x+pad-4)] + 8 * imgPad[(y+pad-size)*W+(x+pad-3)] +
				10 * imgPad[(y+pad-size)*W+(x+pad-2)] + 11 * imgPad[(y+pad-size)*W+(x+pad-1)] + 11 * imgPad[(y+pad-size)*W+(x+pad)] + 11 * imgPad[(y+pad-size)*W+(x+pad+1)] +
				10 * imgPad[(y+pad-size)*W+(x+pad+2)] + 8 * imgPad[(y+pad-size)*W+(x+pad+3)] + 7 * imgPad[(y+pad-size)*W+(x+pad+4)] + 5 * imgPad[(y+pad-size)*W+(x+pad+5)] +
				4 * imgPad[(y+pad-size)*W+(x+pad+6)] + 2 * imgPad[(y+pad-size)*W+(x+pad+7)] + 2 * imgPad[(y+pad-size)*W+(x+pad+8)] + 1 * imgPad[(y+pad-size)*W+(x+pad+9)] +
				1 * imgPad[(y+pad-size)*W+(x+pad+10)] -
				1*imgPad[(y+pad+size)*W+(x+pad-10)] - 1*imgPad[(y+pad+size)*W+(x+pad-9)] - 2*imgPad[(y+pad+size)*W+(x+pad-8)] - 2*imgPad[(y+pad+size)*W+(x+pad-7)] -
				4*imgPad[(y+pad+size)*W+(x+pad-6)] - 5 * imgPad[(y+pad+size)*W+(x+pad-5)] - 7 * imgPad[(y+pad+size)*W+(x+pad-4)] - 8 * imgPad[(y+pad+size)*W+(x+pad-3)] -
				10 * imgPad[(y+pad+size)*W+(x+pad-2)] - 11 * imgPad[(y+pad+size)*W+(x+pad-1)] - 11 * imgPad[(y+pad+size)*W+(x+pad)] - 11 * imgPad[(y+pad+size)*W+(x+pad+1)] -
				10 * imgPad[(y+pad+size)*W+(x+pad+2)] - 8 * imgPad[(y+pad+size)*W+(x+pad+3)] - 7 * imgPad[(y+pad+size)*W+(x+pad+4)] - 5 * imgPad[(y+pad+size)*W+(x+pad+5)] -
				4 * imgPad[(y+pad+size)*W+(x+pad+6)] - 2 * imgPad[(y+pad+size)*W+(x+pad+7)] - 2 * imgPad[(y+pad+size)*W+(x+pad+8)] - 1 * imgPad[(y+pad+size)*W+(x+pad+9)] -
				1 * imgPad[(y+pad+size)*W+(x+pad+10)]);
			imdx[y*W+x] = (long int)(1*imgPad[(y+pad-10)*W+(x+pad-size)] + 1*imgPad[(y+pad-9)*W+(x+pad-size)] + 2*imgPad[(y+pad-8)*W+(x+pad-size)] + 2*imgPad[(y+pad-7)*W+(x+pad-size)] +
				4*imgPad[(y+pad-6)*W+(x+pad-size)] + 5 * imgPad[(y+pad-5)*W+(x+pad-size)] + 7 * imgPad[(y+pad-4)*W+(x+pad-size)] + 8 * imgPad[(y+pad-3)*W+(x+pad-size)] +
				10 * imgPad[(y+pad-2)*W+(x+pad-size)] + 11 * imgPad[(y+pad-1)*W+(x+pad-size)] + 11 * imgPad[(y+pad)*W+(x+pad-size)] + 11 * imgPad[(y+pad+1)*W+(x+pad-size)] +
				10 * imgPad[(y+pad+2)*W+(x+pad-size)] + 8 * imgPad[(y+pad+3)*W+(x+pad-size)] + 7 * imgPad[(y+pad+4)*W+(x+pad-size)] + 5 * imgPad[(y+pad+5)*W+(x+pad-size)] +
				4 * imgPad[(y+pad+6)*W+(x+pad-size)] + 2 * imgPad[(y+pad+7)*W+(x+pad-size)] + 2 * imgPad[(y+pad+8)*W+(x+pad-size)] + 1 * imgPad[(y+pad+9)*W+(x+pad-size)] +
				1 * imgPad[(y+pad+10)*W+(x+pad-size)] -
				1*imgPad[(y+pad-10)*W+(x+pad+size)] - 1*imgPad[(y+pad-9)*W+(x+pad+size)] - 2*imgPad[(y+pad-8)*W+(x+pad+size)] - 2*imgPad[(y+pad-7)*W+(x+pad+size)] -
				4*imgPad[(y+pad-6)*W+(x+pad+size)] - 5 * imgPad[(y+pad-5)*W+(x+pad+size)] - 7 * imgPad[(y+pad-4)*W+(x+pad+size)] - 8 * imgPad[(y+pad-3)*W+(x+pad+size)] -
				10 * imgPad[(y+pad-2)*W+(x+pad+size)] - 11 * imgPad[(y+pad-1)*W+(x+pad+size)] - 11 * imgPad[(y+pad)*W+(x+pad+size)] - 11 * imgPad[(y+pad+1)*W+(x+pad+size)] -
				10 * imgPad[(y+pad+2)*W+(x+pad+size)] - 8 * imgPad[(y+pad+3)*W+(x+pad+size)] - 7 * imgPad[(y+pad+4)*W+(x+pad+size)] - 5 * imgPad[(y+pad+5)*W+(x+pad+size)] -
				4 * imgPad[(y+pad+6)*W+(x+pad+size)] - 2 * imgPad[(y+pad+7)*W+(x+pad+size)] - 2 * imgPad[(y+pad+8)*W+(x+pad+size)] - 1 * imgPad[(y+pad+9)*W+(x+pad+size)] -
				1 * imgPad[(y+pad+10)*W+(x+pad+size)]);
			val = (double)(imdx[y*W+x] * imdx[y*W+x] + imdy[y*W+x]*imdy[y*W+x]);
			val = sqrt(val);
			dest[y*W+x] = (long)val;
			if (dest[y*W+x] > max)
				max = dest[y*W+x];
		}
	}
	printf("maxmax %ld\n", max);

	if (flags&FLT_SOBEL_SELECTION)
	{
		max = 0;
		imgHist = (long*)malloc(H*W*sizeof(long));
		memcpy(imgHist,dest,H*W*sizeof(long));
		qsort(imgHist,H*W,sizeof(long),sort_long);
		lowMag = imgHist[3*H*W/4+1];
		highMag = imgHist[99*H*W/100+1];
		printf("low %ld high %ld \n", lowMag, highMag);
		for (cPt = 0; cPt < W*H; cPt++)
		{
			//fprintf(log, "%d ", imgHist[cPt]);
			if (dest[cPt] < lowMag || dest[cPt] > highMag)
			{
				dest[cPt] = 0;
				imdx[cPt] = 0;
				imdy[cPt] = 0;
			}
			if (dest[cPt] > max)
				max = dest[cPt];
		}
		//fprintf(log, "\n\n");
		free(imgHist);
	}

		for (y = 0; y < H; y++)
		{
			for(x = 0; x < W; x++)
			{
				imgSobel[y * W + x] = (unsigned char)((double)dest[y*W+x]/max * 255);
				fprintf(log, "%u ", imgSobel[y*W+x]);
			}
			fprintf(log, "\n");
		}
	if (flags&FLT_SOBEL_SAVEIMAGE)
		SaveBmp8(name,"_bigsobel.bmp", W, H, imgSobel, 4);
	free(imgPad);
	fclose(log);
	return 0;
}

int FLT_Sobel3x3(int* dest,//destination image
				 unsigned char* imgSobel,//gradient magnitude image
				 const unsigned char *img,//original image
				 char *name,//original image name
				 int H,//height 
				 int W,//width
				 int* mask,	//Sobel mask using for gradient calculation
				 int* imdx,	//gradient values in x direction
				 int* imdy,	//gradient values in y direction
				 int flags)
{
	FILE* log = fopen("../data/res/0629/sobel_log.txt", "wb");
	int pad,size,x,y, cPt, max, lowMag, highMag;
	double val;
	int* imgPad, *imgHist;
	memset(dest, 0, H*W*sizeof(unsigned char));
	memset(imdx, 0, H*W*sizeof(int));
	memset(imdy, 0, H*W*sizeof(int));
	memset(imgSobel,(unsigned char)0, H*W*sizeof(unsigned  char));
	//Padding

	pad = 1;
	imgPad = (unsigned int*)malloc((H+2*pad)*(W+2*pad)*sizeof(int));
	
	for (y = 0; y < H; y++)
	{
		for (x = 0; x < W; x++)
		{
			imgPad[(y+pad) * W + (x+pad)] = (int)img[y * W + x];
		}
	}
	for (y = 0; y < pad; y++)
	{
		for (x = pad; x < W + pad; x++)
		{
			imgPad[y * W + x] = (int)img[(2 * pad - y) * W + x];
			//
			imgPad[(H+y) * W + x] = (int)img[(H-y) * W + x];
			//
		}
	}
	for (x = 0; x < pad; x++)
	{
		for (y = pad; y < H+pad; y++)
		{
			imgPad[y * W + x] = (int)img[y * W + (2 * pad-x)];
			imgPad[y * W + W + x] = (int)img[y*W + W - x];
		}
	}
	for (y = 0; y < pad; y++)
	{
		for (x = 0; x < pad; x++)
		{
			imgPad[y * W + x] = (imgPad[y * W + (2*pad-x)] + imgPad[(2*pad - y) * W + x]) / 2;
			imgPad[(H+y) * W + x] = (imgPad[(H+y) * W + (2*pad-x)] + imgPad[(H-y) * W + x]) / 2;
			imgPad[y * W + W + x] = (imgPad[y * W + (W-x)] + imgPad[(pad - y) * W + x+W]) / 2;
			imgPad[(H+y) * W + W + x] = (imgPad[(H+y) * W + (W-x)] + imgPad[(H-y) * W + x+W]) / 2;
		}
	}
	//Gradient X and Y calculation
	/*
	for (y = 0; y < 2*pad; y++)
	{
		for (x = 0; x < 2*pad; x++)
			fprintf(log, "%d ", imgPad[y*W+x]);
		fprintf(log, "\n");
	}*/
	size = 1;
	max = 0;
	
	for (y = 0; y < H; y++)
	{
		for (x = 0; x < W; x++)
		{
			imdy[y*W+x] = mask[0] * imgPad[(y+pad-1)*W+(x+pad-1)] + mask[1] * imgPad[(y+pad-1)*W+(x+pad)] + mask[2] * imgPad[(y+pad-1)*W+(x+pad+1)] + 
				mask[3] * imgPad[(y+pad)*W+(x+pad-1)] + mask[4]*imgPad[(y+pad)*W+(x+pad)] + mask[5]*imgPad[(y+pad)*W+(x+pad+1)] + 
				mask[6] * imgPad[(y+pad+1)*W+(x+pad-1)] + mask[7]*imgPad[(y+pad+1)*W+(x+pad)] + mask[8]*imgPad[(y+pad+1)*W+(x+pad+1)];

			imdx[y*W+x] = mask[0] * imgPad[(y+pad-1)*W+(x+pad-1)] + mask[3] * imgPad[(y+pad-1)*W+(x+pad)] + mask[6] * imgPad[(y+pad-1)*W+(x+pad+1)] + 
				mask[1] * imgPad[(y+pad)*W+(x+pad-1)] + mask[4]*imgPad[(y+pad)*W+(x+pad)] + mask[7]*imgPad[(y+pad)*W+(x+pad+1)] + 
				mask[2] * imgPad[(y+pad+1)*W+(x+pad-1)] + mask[5]*imgPad[(y+pad+1)*W+(x+pad)] + mask[8]*imgPad[(y+pad+1)*W+(x+pad+1)];
			val = (double)(imdx[y*W+x] * imdx[y*W+x] + imdy[y*W+x]*imdy[y*W+x]);
			val = sqrt(val);
			dest[y*W+x] = (int)val;
			if (dest[y*W+x] > max)
				max = dest[y*W+x];
		}
	}
	if (flags&FLT_SOBEL_SELECTION)
	{
		max = 0;
		imgHist = (int*)malloc(H*W*sizeof(int));
		memcpy(imgHist,dest,H*W*sizeof(int));
		qsort(imgHist,H*W,sizeof(int),sort_int);
		lowMag = imgHist[3*H*W/4+1];
		highMag = imgHist[99*H*W/100];
		for (cPt = 0; cPt < W*H; cPt++)
		{
			if (dest[cPt] < lowMag || dest[cPt] > highMag)
			{
				dest[cPt] = 0;
				imdx[cPt] = 0;
				imdy[cPt] = 0;
			}
			if (dest[cPt] > max)
				max = dest[cPt];
		}
		free(imgHist);
	}


	//fprintf(log,"max = %d.\n", max);
	for (y = 0; y < H; y++)
	{
		for(x = 0; x < W; x++)
		{
			imgSobel[y * W + x] = (unsigned char)((double)dest[y*W+x]/max * 255);
	//		fprintf(log, "%u ", imgSobel[y*W+x]);
		}
	//	fprintf(log, "\n");
	}
	if (flags&FLT_SOBEL_SAVEIMAGE)
	{
		SaveBmp8(name,"_sobel.bmp", W, H, imgSobel, 4);
	}
	free(imgPad);
	
	fclose(log);
	return 0;
}