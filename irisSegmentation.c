#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "header.h"
#include "imio.h"

int ibd_IrisSegmentation(CInfo* CI,CInfo* CP,CInfo *rCP, char* name, unsigned char* imgInput,  unsigned char *imgOutput, int H, int W, int angle, FILE* log, int flags)
{
	int res, Wseg,Hseg, i, rbeg;
	int *xpList, *ypList, xp_mean, yp_mean, rp_mean; //pupil refinement lists of points (x- and y-coordinates)
	uint8 *imgBl, *tmp;
	if ((res = BDL_Flash_RemoveFlash(imgInput, W, H)) != 0)
		fprintf(log, "Flash removal error\n");
	/*else
	{
		cname = (char*)malloc(256 * sizeof(char));
		strcpy(cname, name);
		strcat(cname,"_noflash.bmp");
		CreateBmp8(cname, W, H, imgInput, 4);
		free(cname);
	}*/
	imgBl = (uint8*)malloc(H*W*sizeof(uint8));
	tmp = (uint8*)malloc(H*W*sizeof(uint8));

	IPL_FILT_HaussBlur5x5(imgBl,imgInput,W,H);
	memcpy(tmp, imgBl,H*W*sizeof(uint8));
	IPL_FILT_HaussBlur5x5(imgBl,tmp,W,H);

	//Apply the gradient pairs method to the image to find the first boundary
	if ((res = ibd_graparBothBorders(CI, CP, name, imgInput,imgBl, H, W, angle, flags)) != 0)
		fprintf(log, "Gradient pair method fail.\n");
	DRAW_CircleInGray(imgOutput, W, H, CP->xc, CP->yc, CP->r, 200);
	rCP->xc = CP->xc;
	rCP->yc = CP->yc;
	rCP->r = CP->r;
	//Draw the boundary circles over the eye image and save the result
	
	//SaveBmp8(name, "_circ.bmp", W, H, imgOutput, 4);

	//Apply the circular shortest path refinement method to the found pupil circle approximation
	rp_mean = CP->r;
	Wseg = 360;
	xpList = (int*)malloc(Wseg * sizeof(int));
	ypList = (int*)malloc(Wseg * sizeof(int));
	if (CP->r > 150)
		rbeg = CP->r - 30;
	else
		rbeg = 20;
	if ((res = IBD_RefinePupil(xpList,ypList,&rp_mean, rCP, rbeg, rCP->r+20, name, imgBl, H, W, flags)) != 0)
	{
		fprintf(log, "Error: Refinement fail\n");
		printf("Error: Refinement fail\n");
	}
	else
	{
		xp_mean = 0;
		yp_mean = 0;
		for (i = 0; i < Wseg; i++)
		{
			xp_mean += xpList[i];
			yp_mean += ypList[i];
//			fprintf(log, "%d ", xpList[i]);
		}
//		fprintf(log, "\n");
//		for (i = 0; i < Wseg; i++)
//			fprintf(log, "%d ", ypList[i]);
		xp_mean /= Wseg;
		yp_mean /= Wseg;
		rCP->xc = xp_mean+1;
		rCP->yc = yp_mean+1;
		rCP->r = rp_mean;
		//Draw the refined pupil boundary over the eye image with the previously saved results
		DRAW_SequenceInGray(imgOutput,W,H,xpList,ypList,Wseg, 255);
	}

	/*if ((res = IBD_RefineIris(rCP, CI, (const unsigned char*)imgInput, name, H, W, (BDL_PUPREF_SAVEPOLAR))) != 0)
	{
	//	fprintf(log, "Error: Iris refinement fail\n");
		printf("Error: Iris refinement fail\n");
	}*/
	
	free(xpList);
	free(ypList);
	free(imgBl);
	free(tmp);
	return 0;
}