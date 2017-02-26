#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "header.h"
#include "imio.h"
#include "IrisSegmentation.h"


int IS_Initialize(SSegmentationResult* pResult)
{
	if (pResult == NULL)
	{
		printf("\nError: Result struct was not initialized.");
		return (-1);
	}
	
	pResult->IrisData = (SIrisData*)malloc(sizeof(SIrisData));
	pResult->PupilData = (SPupilData*)malloc(sizeof(SPupilData));

	pResult->IrisData->sCircle = (SCircleData*)malloc(sizeof(SCircleData));
	pResult->IrisData->sCircle->xc = -1;
	pResult->IrisData->sCircle->yc = -1;
	pResult->IrisData->sCircle->r = -1;
	pResult->IrisData->quality = -1.0;

	
	pResult->PupilData->sCircle = (SCircleData*)malloc(sizeof(SCircleData));
	pResult->PupilData->sCircle->xc = -1;
	pResult->PupilData->sCircle->yc = -1;
	pResult->PupilData->sCircle->r = -1;
	
	pResult->PupilData->sRCircle = (SCircleData*)malloc(sizeof(SCircleData));
	pResult->PupilData->sRCircle->xc = -1;
	pResult->PupilData->sRCircle->yc = -1;
	pResult->PupilData->sRCircle->r = -1;
	pResult->PupilData->sContour = (sPoint*)malloc(CONTOUR_QUANTIZATION * sizeof(sPoint));
	pResult->PupilData->quality = -1.0;

	pResult->NormalizedIris = NULL; // Not defined what to do here yet

	return 0;
}

int IS_Deinitialize(SSegmentationResult* pResult)
{
	if (pResult == NULL || pResult->IrisData == NULL || pResult->PupilData == NULL)
		return (-1);
	
	free(pResult->PupilData->sContour);
	free(pResult->PupilData->sRCircle);
	free(pResult->PupilData->sCircle);
	free(pResult->PupilData);

	free(pResult->IrisData->sCircle);
	free(pResult->IrisData);

	// uncomment on update
	// free(pResult->NormalizedIris);

	return 0;
}



//TODO: log - with #define,
//TODO: change the pipeline
//TODO: modify angle params and flags
//TODO: add edgemap and pass in into segmentation functions to calculate it once and for all
int IrisSegmentation(SSegmentationResult* Result, //Segmentation result
	uint8 *imgInput, int H, int W,
	int angle, int flags)
{
	if (Result == NULL)
	{
		printf("\nError: IrisSegmentation: Segmentation result was not initialized.");
		return -1;
	}

	int res, Wseg, Hseg, i, rbeg;
	int xp_mean, yp_mean, pupRadEstimate; //pupil refinement lists of points (x- and y-coordinates)
	uint8 *imgBl, *tmp;
	/*
	if ((res = BDL_Flash_RemoveFlash(imgInput, W, H)) != 0)
	fprintf(log, "Flash removal error\n");
	else
	{
	cname = (char*)malloc(256 * sizeof(char));
	strcpy(cname, name);
	strcat(cname,"_noflash.bmp");
	CreateBmp8(cname, W, H, imgInput, 4);
	free(cname);
	}*/


	// pre-processing stage
	imgBl = (uint8*)malloc(H*W * sizeof(uint8));
	tmp = (uint8*)malloc(H*W * sizeof(uint8));

	IPL_FILT_HaussBlur5x5(imgBl, imgInput, W, H);
	memcpy(tmp, imgBl, H*W * sizeof(uint8));
	IPL_FILT_HaussBlur5x5(imgBl, tmp, W, H);

	//Apply the gradient pairs method to the image to find the first boundary
	if ((res = IS_ApproxBothBorders(Result, imgInput, H, W, angle, flags)) != 0)
	{
		printf("Error: Gradient pair method fail.\n");
		return -2;
	}
	
	//DRAW_CircleInGray(imgOutput, W, H, CP->xc, CP->yc, CP->r, 200);
	//rCP->xc = CP->xc;
	//rCP->yc = CP->yc;
	//rCP->r = CP->r;

	//Apply the circular shortest path refinement method to the found pupil circle approximation
	if ((res = IS_RefinePupil(Result, imgBl, H, W, flags)) != 0)
	{
		fprintf(stderr, "Error: Refinement fail\n");
		return -1;
	}
		
	/*if ((res = IBD_RefineIris(rCP, CI, (const unsigned char*)imgInput, name, H, W, (BDL_PUPREF_SAVEPOLAR))) != 0)
	{
	//	fprintf(log, "Error: Iris refinement fail\n");
	printf("Error: Iris refinement fail\n");
	}*/
	free(imgBl);
	free(tmp);
	return 0;
}

static int sort_edgepnt(const void* e1, const void* e2)
{
	return ((SEdgePnt*)e1)->dir - ((SEdgePnt*)e2)->dir;
}

static int sort_img8U(const void* e1, const void* e2)
{
	return *((unsigned char*)e1) - *((unsigned char*)e2);
}

int IS_ApproxBothBorders(SSegmentationResult *Result, const uint8 *img, int H, int W, int dPhi, int flags)
{
	//return -1;
	clock_t time;
	int nGp, dir, cGp, gval, val, scalar, x, y, x0, y0, cPt, nBegPairIdx, nEndPairIdx, tint, cPar, D, Dt, rmn, rmx, nC, cC, vtype = 0;
	int dx, dy, dr1, dr2;
	int xt, yt, dx1, dx2, dy1, dy2;
	float t, r1, r2, xc, yc, dr;
	double sigma, tlow, thigh, vrmax;
	int xmax = -1, ymax = -1, vmax, rmax = -1, maxRad;
	SCircleData CI, CP;

	uint8 med, max;
	FILE* rhout;

	SEdgePnt* pnEdgePnts;
	//SCircleCnd* pnCircleCnd;
	
	int *pnAcc, *pnAccBl;
	uint8* pucEdge = NULL;
	uint8 *imgsort, *imgMorph;
	void* canny_buf = NULL;
	double* radHist, *radHistSmooth;
	int16 *pgx_calc, *pgy_calc;

	int canny_sz;
	int sz;
	int res = 0;
	sz = 0;
	sz += H*W * sizeof(pnEdgePnts[0]) +  // edge array
		H*W * sizeof(pnAcc[0]) +     // accum
		H*W * sizeof(pnAccBl[0]);// +     // bluraccum
		//H*W * sizeof(pnCircleCnd[0]); // circle params 


	pnEdgePnts = (SEdgePnt*)malloc(H*W * sizeof(pnEdgePnts[0]));
	pnAcc = (int*)malloc(H*W * sizeof(pnAcc[0]));
	pnAccBl = (int*)malloc(H*W * sizeof(pnAccBl[0]));
	//pnCircleCnd = (SCircleCnd*)malloc(2 * H*W * sizeof(pnCircleCnd[0]));
	
	pucEdge = (uint8*)malloc(H*W * sizeof(pucEdge[0]));
	canny_buf = (void*)malloc(4 * sz);

	//median search for canny tresholding
	imgsort = (uint8*)malloc(H*W * sizeof(imgsort[0]));
	memcpy(imgsort, img, sizeof(img[0])*H*W);
	
	//TODO: Estimate image median fast
	qsort(imgsort, H*W, sizeof(uint8), sort_img8U);


	med = imgsort[H*W / 2 + 1];
	max = imgsort[H*W - 1];
	free(imgsort);

	//Canny edge detection and intensity gradients calculating
	//automated canny tresholds calc
	t = (double)med / (int)max;
	printf("t = %f, max = %u\n", t, max);
	tlow = max(0.55, 0.8 * t);
	thigh = min(.8f, 1.33 * t);
	//tlow = .55f;
	//thigh = .7f;
	sigma = 2.f + 2.5*t;
	/*else
	{
	tlow = .4f;
	thigh = .6f;
	sigma = 7.f;
	}*/

	canny_sz = 4 * sz; //5222400
	if ((res = IPL_FILT_Canny(pucEdge, img, &pgx_calc, &pgy_calc, W, H,
		sigma, //5.f,
		tlow,
		thigh,
		canny_buf, &canny_sz, "name")) != ERR_OK)					//original: 4.f .7f .8f
	{
		printf("Error %d: Canny\n", res);
		//return -1;		
	}

	/*
	if (flags&BDL_PGRAD_SAVECANNY)
	SaveBmp8("D:/data/res/0224/", "_canny.bmp", W, H, pucEdge, 4);
	*/
	SaveBmp8("D:/data/res/0224/", "_canny.bmp", W, H, pucEdge, 4);

	time = clock();
	// Collect gradient point to array
	nGp = 0;
	for (y = 10; y < H - 10; y++)
	{
		for (x = 10; x < W - 10; x++)
		{
			if (pucEdge[cPt = y*W + x])
			{
				dir = (int)(.5 + 180.*atan2((float)pgy_calc[cPt], (float)pgx_calc[cPt]) / PI);
				pnEdgePnts[nGp].dir = dir;
				pnEdgePnts[nGp].x = (int16)x;
				pnEdgePnts[nGp].y = (int16)y;
				pnEdgePnts[nGp].gx = pgx_calc[cPt];
				pnEdgePnts[nGp].gy = pgy_calc[cPt];
				if (dir > -45 && dir < 45 || dir < -135 || dir > 135)
					nGp++;
			}
		}
	}
	printf("nGp: %d\n", nGp);
	time = clock() - time;
	printf("GP extraction: %f sec.\n", (float)time / CLOCKS_PER_SEC);

	//sorting the gradient points array
	qsort(pnEdgePnts, nGp, sizeof(pnEdgePnts[0]), sort_edgepnt);
	printf("Angle values in [%d,%d]\n", pnEdgePnts[0].dir, pnEdgePnts[nGp-1].dir);
	
	//duplicate the array
	memcpy(pnEdgePnts + nGp, pnEdgePnts, sizeof(pnEdgePnts[0])*nGp);

	//clean the accumulator
	for (cGp = 0; cGp<nGp; cGp++)
		pnEdgePnts[cGp + nGp].dir = pnEdgePnts[cGp].dir + 360;

	// clean separated accumulators for center and radius
	memset(pnAcc, 0, W*H * sizeof(pnAcc[0]));

	// prime counters
	tint = pnEdgePnts[0].dir + dPhi - ANGDELT;
	for (nBegPairIdx = 0; ((nBegPairIdx<nGp) &&
		(pnEdgePnts[nBegPairIdx].dir<tint)); nBegPairIdx++);

	tint = pnEdgePnts[0].dir + dPhi + ANGDELT;
	for (nEndPairIdx = nBegPairIdx; ((nEndPairIdx<nGp) &&
		(pnEdgePnts[nEndPairIdx].dir<tint)); nEndPairIdx++);

	nC = 0;
	
	time = clock();
	// main processing
	for (cGp = 0; cGp<nGp; cGp++)
	{
		for (cPar = nBegPairIdx + cGp; cPar < nEndPairIdx + cGp; cPar++)
		{
			D = pnEdgePnts[cPar].gx * pnEdgePnts[cGp].gy - pnEdgePnts[cPar].gy * pnEdgePnts[cGp].gx;
			Dt = pnEdgePnts[cPar].gx * (pnEdgePnts[cPar].y - pnEdgePnts[cGp].y) - pnEdgePnts[cPar].gy * (pnEdgePnts[cPar].x - pnEdgePnts[cGp].x);
			if (D == 0)
				continue;
			
			t = (float)Dt / D;
			dx1 = (int)(pnEdgePnts[cGp].gx * t);
			dy1 = (int)(pnEdgePnts[cGp].gy * t);
			x = pnEdgePnts[cGp].x + dx1;
			y = pnEdgePnts[cGp].y + dy1;
			
			if (x < 1 || y < 1 || x >= W - 10 || y >= H - 10)
				continue;
			
			dx2 = pnEdgePnts[cPar].x - x;
			dy2 = pnEdgePnts[cPar].y - y;
			
			r1 = sqrt(dx1*dx1 + dy1*dy1);
			r2 = sqrt(dx2*dx2 + dy2*dy2);
			
			if (abs(r1 - r2) >= RADDELT || max(r1, r2) < MINRAD)
				continue;
			//Ds = pnEdgePnts[cGp].gx * (pnEdgePnts[cPar].y-pnEdgePnts[cGp].y) - pnEdgePnts[cGp].gy * (pnEdgePnts[cPar].x-pnEdgePnts[cGp].x);
			//s = (double)Ds / (double)D;
			//pnCircleCnd[nC].xc = (int)xt;
			//pnCircleCnd[nC].yc = (int)yt;
			//pnCircleCnd[nC].r = (r1 + r2) / 2;
			//nC++;
			pnAcc[y*W + x]++;
		}
	}
	time = clock() - time;
	printf("Paired Gradient voting: %f sec.\n", (float)time / CLOCKS_PER_SEC);

	maxRad = min(W, H) / 2;
	// blur accumulator
	IPL_FILT_HaussBlur3x3_Int(pnAccBl, pnAcc, W, H, 0);
	//memcpy(pnAcc, pnAccBl, W*H * sizeof(int));
	//IPL_FILT_HaussBlur3x3_Int(pnAccBl, pnAccum, W, H, 0);

	// find maximum
	vmax = xmax = ymax = 0;
	for (y = 1; y<H - 1; y++)
		for (x = 1; x<W - 1; x++)
			if (vmax < pnAccBl[y*W + x])
				vmax = pnAccBl[(ymax = y)*W + (xmax = x)];
	CI.xc = xmax + 1;
	CI.yc = ymax + 1;
	if (vmax > 0)
	{
		uint8* pnAccOut = (uint8*)malloc(H * W * sizeof(uint8));
		// Cleaning borders for accumulator
		for (x = 0; x < max(W, H); x++)
		{
			pnAccOut[x] = 0;
			pnAccOut[(H - 1)*W + x] = 0;
			if (x < min(W, H))
			{
				pnAccOut[x*(W - 1)] = 0;
				pnAccOut[x*(W - 1) + H - 1] = 0;
			}
		}
		// Intensity normalization
		for (y = 1; y < H - 1; y++)
			for (x = 1; x < W - 1; x++)
				pnAccOut[y*W + x] = (uint8)((double)pnAccBl[y*W + x] / vmax * 255);
		/*
		if (flags&BDL_PGRAD_SAVEACC)
		{
			SaveBmp8(name, "_acc1.bmp", W, H, pnAccOut, 4);
		}
		*/
		free(pnAccOut);
	}

	vrmax = 0;
	vmax = rmax = 0;
	if (flags&BDL_PGRAD_CALCRADIUS)
	{
		rmx = max(H, W) / 2;
		radHist = (double*)malloc((H + W) * sizeof(double));
		radHistSmooth = (double*)malloc((H + W) * sizeof(double));
		for (x = 0; x < 400; x++) radHist[x] = 0.0;

		for (cGp = 0; cGp<nGp; cGp++)
		{
			dx1 = (pnEdgePnts[cGp].x - xmax);
			dy1 = (pnEdgePnts[cGp].y - ymax);
			xt = pnEdgePnts[cGp].gx;
			yt = pnEdgePnts[cGp].gy;

			//TODO: Insert fast inverse sqrt here
			t = (dx1*pnEdgePnts[cGp].gx + dy1 * pnEdgePnts[cGp].gy) / sqrt((float)(xt * xt + yt * yt));
			r1 = sqrt(dx1*dx1 + dy1*dy1);
			t = abs(t) / r1;
			if (r1 >= MINRAD && r1 < rmx && t > 0.96)
				radHist[(int)(r1 + .5)]++;
		}

		// Histogram normalization
		for (x = 1; x < rmx; x++) radHist[x] = radHist[x] / x;
		//Histogram Haussian filtering (twice for accuracy)
		IPL_HIST_Blur_double(radHistSmooth, (const double*)radHist, rmx, 5);
		for (x = 1; x < rmx; x++) radHist[x] = radHistSmooth[x];
		IPL_HIST_Blur_double(radHistSmooth, (const double*)radHist, rmx, 5);

		for (x = MINRAD; x<rmx; x++)
			if (vrmax < radHistSmooth[x])
				vrmax = radHistSmooth[rmax = x];
		CI.r = rmax;
	}

	/*
	// SECOND BORDER LOCALIZATION

	//clean accumulator
	memset(pnAcc, 0, W*H * sizeof(pnAcc[0]));

	dr2 = CI.r / 5;
	dr2 = dr2*dr2;
	//CI->yc = H - CI->yc;
	for (cC = 0; cC < nC; cC++)
	{
		dx = pnCircleCnd[cC].xc - CI->xc;
		dy = pnCircleCnd[cC].yc - CI->yc;
		dr1 = dx*dx + dy*dy;
		if (dr1 <= dr2 && (pnCircleCnd[cC].r < 3 * CI->r / 4 || pnCircleCnd[cC].r > 4 * CI->r / 3))
		{
			pnAccum[pnCircleCnd[cC].yc*W + pnCircleCnd[cC].xc]++;
			//radHist[pnCircleCnd[cC].r] *= pnCircleCnd[cC].r;
			//radHistSmooth[pnCircleCnd[cC].r]++;
			//radHist[pnCircleCnd[cC].r] /= pnCircleCnd[cC].r;
			//printf("%d \n", pnCircleCnd[cC].r);
		}
	}
	// blur accumulator
	IPL_FILT_HaussBlur3x3_Int(pnAccBl, pnAcc, W, H, 0);
	memcpy(pnAccum, pnAccBl, H*W * sizeof(int));
	IPL_FILT_HaussBlur3x3_Int(pnAccBl, pnAccum, W, H, 0);

	// find maximum
	vmax = xmax = ymax = 0;
	for (y = 10; y < H - 10; y++)
	{
		for (x = 10; x < W - 10; x++)
		{
			if (vmax < pnAccBl[y*W + x])
				vmax = pnAccBl[(ymax = y)*W + (xmax = x)];
		}
	}
	//set treshold for decentration
	dr1 = CI->r / 5;
	dr1 = dr1 * dr1;

	if ((ymax - CI->yc)*(ymax - CI->yc) + (xmax - CI->xc)*(xmax - CI->xc) < dr1)
	{
		CP->xc = xmax;
		CP->yc = ymax;
	}
	else
	{
		CP->xc = CI->xc;
		CP->yc = CI->yc;
	}

	if (flags&BDL_PGRAD_CALCRADIUS)
	{

		for (dr1 = 0; dr1 <= CI->r / 7; dr1++)
			radHistSmooth[dr1] = 0.0;
		for (dr1 = 0; dr1 <= MINRAD; dr1++)
			radHistSmooth[dr1] = 0.0;
		for (dr1 = 3 * CI->r / 4; dr1 < 4 * CI->r / 3; dr1++)
			radHistSmooth[dr1] = 0.0;

		rmx = max(H, W) / 2;

		for (dr1 = 0; dr1 < rmx; dr1++)
			radHist[dr1] = radHistSmooth[dr1];
		IPL_HIST_Blur_double(radHistSmooth, (const double*)radHist, rmx, 3);
		//'Haussian' filtering for radius search
		vrmax = 0.0;
		for (x = rmx - 1; x>MINRAD; x--)
			if (vrmax < radHistSmooth[x])
				vrmax = radHistSmooth[rmax = x];
		CP->r = rmax;
		
	}
	*/

	Result->IrisData->sCircle->r = CI.r;
	Result->IrisData->sCircle->xc = CI.xc;
	Result->IrisData->sCircle->yc = CI.yc;
	free(radHist);
	free(radHistSmooth);

	free(pnEdgePnts);
	free(pnAcc);
	free(pnAccBl);
	//free(pnCircleCnd);
	free(pucEdge);
	free(canny_buf);
	return 0;
}

int IS_RefinePupil(SSegmentationResult *sResult, uint8 *img, int H, int W, int flags)
{
	char *cname;
	uint8 *imgEdge, *imgPolar, *imgPolarBl, *imgCrop, *imgGrad;
	void *canny_buf = NULL;
	int Wseg, Hseg, i, rad, MaxRad, MinRad, angle, xbeg, ybeg, xend, yend, Hc, Wc, x, y, xc, yc, med, canny_sz, res;
	int gx, gy, ax, ay, *radList, max;
	int *dest, *imdx, *imdy;
	int mask[9] = { 3,10,3,0,0,0,-3,-10,-3 };

	double t, tlow, thigh, sigma, val;
	int16 *pgx_calc, *pgy_calc;
	//FILE *fout = fopen("../data/res/edgeMapVals.txt", "w");

	//Cropping the pupil for better performance
	xbeg = max(sResult->PupilData->sCircle->xc - sResult->PupilData->sCircle->r - W / 10, 1);
	ybeg = max(sResult->PupilData->sCircle->yc - sResult->PupilData->sCircle->r - H / 10, 1);

	xend = min(sResult->PupilData->sCircle->xc + sResult->PupilData->sCircle->r + W / 10, W - 1);
	yend = min(sResult->PupilData->sCircle->yc + sResult->PupilData->sCircle->r + H / 10, W - 1);

	if (yend < 0 || xbeg >= W || ybeg >= H || xend < 0)
	{
		fprintf(stderr, "\nError: RefinePupil: Iris crop wrong rectangle.\n");
		return -1;
	}
	
	Hc = yend - ybeg + 1;
	Wc = xend - xbeg + 1;
	
	xc = sResult->PupilData->sCircle->xc - xbeg;
	yc = sResult->PupilData->sCircle->yc - ybeg;

	imgCrop = (uint8*)malloc(Hc * Wc * sizeof(uint8));
	imgGrad = (uint8*)malloc(Hc * Wc * sizeof(uint8));
	imgEdge = (uint8*)malloc(Hc * Wc * sizeof(uint8));
	dest = (int*)malloc(Hc * Wc * sizeof(int));
	imdx = (int*)malloc(Hc * Wc * sizeof(int));
	imdy = (int*)malloc(Hc * Wc * sizeof(int));

	memset(dest, 0, Hc*Wc * sizeof(int));
	memset(imgCrop, 0, Hc*Wc);
	memset(imgGrad, 0, Hc*Wc);
	canny_buf = (void*)malloc(20 * H * W * sizeof(unsigned char));

	for (y = ybeg; y <= yend; y++)
	{
		for (x = xbeg; x <= xend; x++)
		{
			imgCrop[(y - ybeg)* Wc + (x - xbeg)] = img[y * W + x];
		}
	}

	med = IPL_HIST_mdn_pixE(imgCrop, Hc*Wc);
	//printf("med_crop=%d\n", med);
	t = (double)med / 255;
	//printf("t = %f\n", t);
	tlow = 0.8 * t;
	thigh = 1.33 * t;
	//tlow = .55f;
	//thigh = .7f;
	sigma = 2.f + 3 * t;
	canny_sz = 20 * H * W;//1049665; //Костыль

						  //FLT_Sobel3x3(dest, imgEdge, imgCrop, name, Hc, Wc, mask, imdx, imdy,(FLT_SOBEL_SELECTION));


	if ((res = IPL_FILT_Canny(imgEdge, imgCrop, &pgx_calc, &pgy_calc, Wc, Hc, sigma, tlow, thigh,
		canny_buf, &canny_sz, "_")) != ERR_OK)					//original: 4.f .7f .8f
	{
		printf("Error %d: Canny, refinement stage fault\n", res);
		free(imgCrop);
		free(imgGrad);
		free(imgEdge);
		free(canny_buf);
		return -1;
	}

	max = 0;
	for (y = 1; y < Hc - 1; y++)
	{
		for (x = 1; x < Wc - 1; x++)
		{

			i = Wc * y + x;
			ax = x - xc;
			ay = y - yc;
			gx = (int)pgx_calc[i];
			gy = (int)pgy_calc[i];
			val = (double)(ax * gx + ay * gy);
			val /= (int)(sqrt((double)(ax*ax + ay*ay)));
			val /= (int)(sqrt((double)(gx*gx + gy*gy)));
			//(int)(.5+180.*atan2((double)pgy_calc[cPt],(double)pgx_calc[cPt])/PI)
			//imgGrad[i] = (int)((double)sqrt(gx*gx + gy*gy));
			t = 180. * acos(val) / PI;
			if (t > 10)
			{
				imgEdge[i] = 0;
				//printf("%d %d\n", gx,gy);
			}
			//fprintf(fout, "%d ", imgGrad[i]);
		}
		//fprintf(fout, "\n");
	}
	//fclose(fout);

	/*
	if (flags&BDL_PUPREF_SAVECROP)
	SaveBmp8(name, "_cropcanny.bmp", Wc, Hc, imgEdge, 4);
	*/
	if (sResult->PupilData->sCircle->r > 150)
		MinRad = sResult->PupilData->sCircle->r - 20;
	else
		MinRad = 20;
	MaxRad = sResult->PupilData->sCircle->r + 20;

	Wseg = 360;
	Hseg = MaxRad - MinRad + 1;
	radList = (int*)malloc((Wseg + 1) * sizeof(int));
	printf("improving\n");

	imgPolar = (uint8*)malloc(Hseg*Wseg * sizeof(uint8));
	memset(imgPolar, 0, Hseg*Wseg * sizeof(uint8));


	for (rad = MinRad; rad <= MaxRad; rad++)
	{
		for (angle = 0; angle < Wseg; angle++)
		{
			i = (int)((double)rad * sin((double)(angle) / 180 * PI)) * Wc + (int)((double)rad * cos((double)(angle) / 180 * PI));
			imgPolar[(rad - MinRad)*Wseg + angle] = 255 - imgEdge[Wc * yc + xc + i];
		}
	}
	/*
	if (flags&BDL_PUPREF_SAVEPOLAR)
	SaveBmp8(name, "_polar_refinemap.bmp", Wseg, Hseg, imgPolar, 4);
	*/
	for (rad = 0; rad < Hseg; rad++)
	{
		for (angle = 0; angle < Wseg; angle++)
			if (imgPolar[rad * Wseg + angle] > 0)
				imgPolar[rad * Wseg + angle] = 1;
	}
	if ((res = pupilCircSP(radList, imgPolar, Hseg, Wseg, sResult->PupilData->sCircle->r - MinRad)) != 0)
	{
		free(imgPolar);
		free(radList);
		free(imgCrop);
		free(imgGrad);
		free(imgEdge);
		free(canny_buf);
		printf("Error %d: Circular shortest path obtaining failed\n", res);
		return -1;
	}
	sResult->PupilData->sRCircle->r = 0;
	sResult->PupilData->sRCircle->xc = 0;
	sResult->PupilData->sRCircle->yc = 0;
	x = y = 0;
	for (i = 0; i < Wseg; i++)
	{
		radList[i] += MinRad;
		sResult->PupilData->sRCircle->r += radList[i];
		//printf("%d ", radList[i]);
		sResult->PupilData->sContour[i].x = sResult->PupilData->sCircle->xc + (int)((double)radList[i] * cos((double)(i) / 180 * PI));
		sResult->PupilData->sRCircle->xc += sResult->PupilData->sContour[i].x;
		sResult->PupilData->sContour[i].y = sResult->PupilData->sCircle->yc + (int)((double)radList[i] * sin((double)(i) / 180 * PI));
		sResult->PupilData->sRCircle->yc += sResult->PupilData->sContour[i].y;
		//ypList[i] = CP->yc + (int)((double)radList[i] * sin((double)(i) / 180 * PI));
	}
	sResult->PupilData->sRCircle->r /= Wseg;
	sResult->PupilData->sRCircle->xc /= Wseg;
	sResult->PupilData->sRCircle->yc /= Wseg;

	free(radList);
	free(imgPolar);
	free(canny_buf);
	free(imgCrop);
	free(imgGrad);
	free(imgEdge);
	free(dest);
	free(imdx);
	free(imdy);
	return 0;
}