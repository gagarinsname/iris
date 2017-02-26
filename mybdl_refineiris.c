#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "header.h"
#include "imio.h"
#include "ipl.h"

static int sort_int(const void* e1, const void* e2)
{
	return *((int*)e1) - *((int*)e2);
}


int IBD_RefineIris(CInfo *CP, CInfo *CI, const unsigned char *img, char* name, int H, int W, int flags)
{
	/*
	FILE* rh = fopen("radhist.txt", "wb");
	long res, x,y,xc, yc, gx, gy, val, px, py, rmax, rsize;
	double scalar, rI, rval, ang;
	uint8 *imgSobel, *imgMorph;
	long *imdx, *imdy, *imGrad,*gHist;
	double *radHistSmooth, *radHist;

	rsize = min(H,W);
	rsize /= 2;

	imgMorph = (unsigned char*)malloc(H*W*sizeof(unsigned char));
	Dilate3x3Cross(imgMorph, img, W, H);
	memcpy(img, imgMorph,H*W*sizeof(unsigned char));
	memset(imgMorph,0,H*W*sizeof(unsigned char));
	IPL_FILT_HaussBlur5x5(imgMorph,img,W,H);
	SaveBmp8(name,"_morph.bmp", W, H, imgMorph, 4);
	
	imdx = (long*)malloc(H*W*sizeof(long));
	imdy = (long*)malloc(H*W*sizeof(long));
	imGrad = (long*)malloc(H*W*sizeof(long));
	imgSobel = (unsigned char*)malloc(H*W*sizeof(unsigned char));
	radHist = (double*)malloc((rsize+2)*sizeof(double));
	radHistSmooth = (double*)malloc(rsize*sizeof(double));
	for (x = 0; x < rsize; x++) radHist[x] = 0.0;
	if ((res = FLT_Sobel10x10(imGrad,imgSobel, img, name, H, W, imdx, imdy, 0)) == 0)
	{
		xc = CI->xc;
		yc = CI->yc;
		printf("ctr: %d %d \n", xc,yc);
		for (y = 10; y < H-10; y++)
		{
			for(x = 10; x < W-10; x++)
			{
				px = x - xc;
				py = y - yc;
				if (px*px + py*py > 16 * CP->r*CP->r / 9 && (imdx[y*W+x] != 0 || imdy[y*W+x] != 0))
				{
					gx = imdx[y*W+x];
					gy = imdy[y*W+x];
					scalar = (double)(px*gx + py*gy) / sqrt((double)(gx*gx + gy*gy));
					rI = sqrt((double)(px*px + py*py));
					scalar /= rI;
					ang = 180. * acos(scalar)/PI;

					if (abs(ang) > 30)
						imgSobel[y*W+x] = 0;
					
					if (rI < rsize)
						radHist[(int)(rI+.5)]+=1.0;
					else
						imgSobel[y*W+x] = 0;
				}
				else
					imgSobel[y*W+x] = (unsigned char)0;
			}
		}
		SaveBmp8(name,"_bigsobel.bmp", W, H, imgSobel, 4);
		
		rval = 0.0;
		rmax = CI->r;
		for (x = 1; x < rsize; x++) radHist[x] = radHist[x]/x;
		IPL_HIST_Blur_double(radHistSmooth, (const double*) radHist, rsize, 5);
		for (x = 1; x < rsize; x++) radHist[x] = radHistSmooth[x];
		IPL_HIST_Blur_double(radHistSmooth, (const double*) radHist, rsize, 5);
		for(x = 0; x < rsize; x++)
		{
			if (radHistSmooth[x] > rval)
			{
				rval = radHistSmooth[x];
				rmax = x;
			}
			fprintf(rh, "%1.2f ", radHistSmooth[x]);
		}
		printf("iris radius: %d\n", rmax);
		//CI->r = rmax;
	}
	else
		return -1;
	free(imgMorph);
	free(imdx);
	free(imdy);
	free(imGrad);
	free(imgSobel);
	free(radHist);
	free(radHistSmooth);
	fclose(rh);
	*/
	return 0;
}