#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "header.h"
#include "imio.h"
#include "IrisSegmentation.h"

double frac(double val) {
	return val - floor(val);
}

static int sort_edgepnt(const void* e1, const void* e2)
{
	return ((SEdgePnt*)e1)->dir - ((SEdgePnt*)e2)->dir;
}

static int sort_img8U(const void* e1, const void* e2)
{
	return *((unsigned char*)e1) - *((unsigned char*)e2);
}

static int sort_int(const void *e1, const void *e2)
{
	return *((int*)e1) - *((int*)e2);
}


int IS_Initialize(SSegmentationResult* pResult, char* name)
{
	if (pResult == NULL)
	{
		printf("\nError: Result struct was not initialized.");
		return (-1);
	}
	
	pResult->name = (char*)malloc(512 * sizeof(char));
	strcpy(pResult->name, name);

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
	
	free(pResult->name);

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

	// pre-processing stage
	imgBl = (uint8*)malloc(H*W * sizeof(uint8));
	tmp = (uint8*)malloc(H*W * sizeof(uint8));

	Dilate3x3Cross(imgBl, imgInput, W, H);
	memcpy(tmp, imgBl, H*W * sizeof(uint8));
	Dilate3x3Cross(imgBl, tmp, W, H);
	//SaveBmp8("D:/data/res/0615/1.bmp", "cross", W, H, imgBl, 4);
	memcpy(tmp, imgBl, H*W * sizeof(uint8));
	IPL_FILT_HaussBlur5x5(imgBl, tmp, W, H);
	

	//Apply the gradient pairs method to the image to find the first boundary
	if ((res = IS_ApproxBothBorders(Result, imgInput, H, W, angle, flags)) != 0)
	{
		fprintf(stderr, "Error: Gradient pair method fail.\n");
		return -2;
	}

	
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

/* ������� ���������� ����� ����������� ��� ������������ ����������� image � ����� ������ �������� size */
uint8 otsuThreshold(const unsigned char *image, const int size)
{
	uint8 min = image[0], max = image[0];
	uint8 tmp;
	int temp, temp1;
	int *hist;
	int histSize;
	uint8 *imgSave = (uint8*)malloc(size * sizeof(uint8));

	int alpha, beta, threshold = 0;
	double sigma, maxSigma = -1;
	double w1, a;

	/**** ���������� ����������� ****/
	/* ������ ���������� � ���������� ������� */
	for (int i = 1; i < size; ++i)
	{
		tmp = image[i];
		if (tmp < min)   min = tmp;
		if (tmp > max)   max = tmp;
	}

	histSize = (int)(max - min) + 1;
	if ((hist = (int*) malloc(sizeof(int) * histSize)) == NULL) return -1;

	memset(hist, 0, histSize * sizeof(int));

	/* ������� ������� ����� ��������� */
	for (int i = 0; i < size; ++i)
		++hist[(int)(image[i] - min)];

	/**** ����������� ��������� ****/

	temp = temp1 = 0;
	alpha = beta = 0;
	/* ��� ������� ��������������� �������� ������� ������ */
	for (int i = 0; i <= (int)(max - min); ++i)
	{
		temp += i*hist[i];
		temp1 += hist[i];
	}

	/* �������� ���� ������ ������
	����������� �� ���� ��������� ��� ������ ������, ��� ������� ��������������� ��������� ���������� */
	for (int i = 0; i<(int)(max - min); ++i)
	{
		alpha += i*hist[i];
		beta += hist[i];

		w1 = (double)beta / temp1;
		a = (double)alpha / beta - (double)(temp - alpha) / (temp1 - beta);
		sigma = w1*(1 - w1)*a*a;

		if (sigma>maxSigma)
		{
			maxSigma = sigma;
			threshold = i;
		}
	}
	free(hist);
	
	return (uint8)threshold + min;
}

/* Function is used to search for eye center given the filled accumulator array*/
int pgDetectCenter(SCircleData *CI, SCircleData *CP, int* pnAcc, int H, int W, int flags)
{
	int x,y, vmax,xmax,ymax;
	int* pnAccBl = (int*)malloc(H*W * sizeof(int));

	// blur accumulator
	IPL_FILT_HaussBlur3x3_Int(pnAccBl, pnAcc, W, H, 1);

	// find maximum
	vmax = xmax = ymax = 0;
	for (y = 1; y < H - 1; ++y)
		for (x = 1; x < W - 1; ++x)
			if (vmax < pnAccBl[y*W + x])
				vmax = pnAccBl[(ymax = y)*W + (xmax = x)];
	
	printf("Max vote: %d\n", vmax);

	CP->xc = xmax;
	CP->yc = ymax;
	
	CI->xc = xmax;
	CI->yc = ymax;


	if (flags&BDL_PGRAD_SAVEACC)
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
		SaveBmp8(name, "_acc1.bmp", W, H, pnAccOut, 4);
		*/
		free(pnAccOut);
	}
	
	free(pnAccBl);

	if (vmax == 0)
		return -1;

	return 0;
}


int pgEstimateRadii(SCircleData *CP, SCircleData *CI, SEdgePnt* pnEdgePnts, int nGp, int H, int W, int flags)
{
	int dx1, dy1, gx, gy, r, rMax, vrmax, vmax, cGp;
	float t, dr;
	int *radHist, *radHistSmooth;

	rMax = min(min(CP->xc, W - CP->xc), min(CP->yc, H - CP->yc));
	radHist = (int*)malloc((H + W) * sizeof(int));
	radHistSmooth = (int*)malloc((H + W) * sizeof(int));
	memset(radHist, 0, (H + W) * sizeof(int));
	float threshDot = 0.96f * 0.96f;
	int minRad2 = MINRAD * MINRAD, maxRad2 = rMax*rMax;

	for (cGp = 0; cGp<nGp; ++cGp)
	{
		dx1 = pnEdgePnts[cGp].x - CI->xc;
		dy1 = pnEdgePnts[cGp].y - CI->yc;
		gx = pnEdgePnts[cGp].gx;
		gy = pnEdgePnts[cGp].gy;

		int dotProduct = dx1*pnEdgePnts[cGp].gx + dy1 * pnEdgePnts[cGp].gy;
		/// sqrtf((float)(gx * gx + gy * gy));
		int dr = dx1*dx1 + dy1*dy1;
		// t = abs(t) / dr;
		if (dr >= minRad2 && dr < maxRad2 && dotProduct * dotProduct > threshDot * dr)
			++radHist[(int)(sqrtf((float)(dr + .5)))];
	}

	// Histogram normalization
	const int scaleConstant = 4096;
	for (r = 1; r < rMax; r++) radHist[r] = scaleConstant * radHist[r] / r;
	//Histogram Haussian filtering (twice for accuracy)
	int blur_hw = 5;
	IPL_HIST_Blur(radHistSmooth, (const int*)radHist, rMax + 1, blur_hw);
	free(radHist);

	vrmax = 0;
	for (r = MINRAD; r < rMax; ++r)
		if (vrmax < radHistSmooth[r])
			vrmax = radHistSmooth[vmax = r];

	int vrmax2 = 0;
	int hRmax = vmax / 7;

	for (r = 0; r < hRmax; ++r)
		radHistSmooth[r] = 0.0;

	hRmax = 4 * vmax / 3;
	int lRmax = 3 * vmax / 4;
	for (r = lRmax; r < hRmax; ++r)
		radHistSmooth[r] = 0.0;

	int vmax2 = 0;
	for (r = MINRAD; r < rMax; ++r)
		if (vrmax2 < radHistSmooth[r])
			vrmax2 = radHistSmooth[vmax2 = r];

	if (vmax > vmax2)
	{
		int sw = vmax;
		vmax = vmax2;
		vmax2 = sw;
	}

	free(radHistSmooth);

	if (vmax > 0)
	{
		CP->r = vmax;
		if (vmax2 > 0)
			CI->r = vmax2;
	}
	else
		return -1;

	return 0;
}


int IS_ApproxBothBorders(SSegmentationResult *Result, const uint8 *img, int H, int W, int dPhi, int flags)
{
	clock_t time;
	int nGp, dir, cGp, gval, val, x, y, x0, y0, cPt, nBegPairIdx, nEndPairIdx, tint, cPar, D, Dt, rmn, rmx, nC, cC, vtype = 0;
	int dx, dy, dr1, dr2, r;
	int xmax = -1, ymax = -1, vmax, rmax = -1, maxRad;
	int xt, yt, dx1, dx2, dy1, dy2;
	float t, r1, r2, xc, yc, dr;
	double sigma, tlow, thigh, tOtsu;

	SCircleData CI, CP;

	uint8 med, max;
	FILE* rhout;

	SEdgePnt* pnEdgePnts;

	int *pnAcc;
	uint8* pucEdge = NULL;
	uint8 *imgsort, *imgMorph;
	void* canny_buf = NULL;
	int16 *pgx_calc, *pgy_calc;

	int canny_sz;
	int sz;
	int res = 0;
	sz = 0;
	sz += H*W * sizeof(pnEdgePnts[0]) +  // edge array
		H*W * sizeof(pnAcc[0]);     // accum
		//H*W * sizeof(pnCircleCnd[0]); // circle params 


	pnEdgePnts = (SEdgePnt*)malloc(H*W * sizeof(pnEdgePnts[0]));
	pnAcc = (int*)malloc(H*W * sizeof(pnAcc[0]));
	// pnAccBl = (int*)malloc(H*W * sizeof(pnAccBl[0]));

	pucEdge = (uint8*)malloc(H*W * sizeof(pucEdge[0]));
	canny_buf = (void*)malloc(4 * sz);

	//median search for canny tresholding
	imgsort = (uint8*)malloc(H*W * sizeof(imgsort[0]));
	memcpy(imgsort, img, sizeof(img[0])*H*W);
	//TODO: Estimate image median fast
	qsort(imgsort, H*W, sizeof(uint8), sort_img8U);

	max = imgsort[H*W - 1];
	free(imgsort);
	tOtsu = otsuThreshold(img, H*W);

	uint8 *imgOtsu = (uint8*)malloc(H*W * sizeof(uint8));
	for (int i = 0; i < H * W; ++i)
		imgOtsu[i] = img[i] > tOtsu ? 0 : 255;
	int *xProj = (int*)malloc(W * sizeof(int));
	int *yProj = (int*)malloc(H * sizeof(int));

	int xBeg = 0, yBeg = 0, xEnd = W, yEnd = H;
	for (int x = 0; x < W; ++x)
	{
		for (int y = 0; y < H; ++y)
		{
			xProj[x] += imgOtsu[x * H + y];
			yProj[y] += imgOtsu[x * H + y];
		}
	} 
	for (int x = 0; x < W; ++x)
	{
		xProj[x] /= H;
		if (xProj != 255)
			break;
		xBeg = x;
	}
	for (int y = 0; y < H; ++y)
		yProj[y] /= W;


	// SaveBmp8(Result->name, "_otsu", W, H, imgOtsu, 4);
	free(imgOtsu);

	printf("[PG] t_otsu = %f\n", (double)tOtsu / 255);

	//Canny edge detection and intensity gradients calculating
	//automated canny tresholds calc
	//Otsu thresholding for Canny edge detection
	tOtsu = (double)(otsuThreshold((const unsigned char*)img, (const int)H*W)) / 255;
	// t = (double)med / (int)max;
	thigh = min(.7f, 1.2 * tOtsu / 640 * W);
	tlow = thigh / 2;
	//tlow = .55f;
	//thigh = .7f;
	sigma = 4.f;

	canny_sz = 4 * sz; //5222400
	if ((res = IPL_FILT_Canny(pucEdge, img, &pgx_calc, &pgy_calc, W, H,
		sigma, //5.f,
		tlow,
		thigh,
		canny_buf, &canny_sz, "name")) != ERR_OK)					//original: 4.f .7f .8f
	{
		printf("[ERROR]: Canny failed %d.\n", res);
		return -1;		
	}
	// SaveBmp8("D:/data/res/0224/file.bmp", "_canny.bmp", W, H, pucEdge, 4);

	time = clock();
	// Collect gradient point to array
	nGp = 0;
	for (y = 10; y < H - 10; ++y)
	{
		for (x = 10; x < W - 10; ++x)
		{
			if (pucEdge[cPt = y*W + x])
			{
				dir = (int)(.5 + 180.*atan2((float)pgy_calc[cPt], (float)pgx_calc[cPt]) / PI);
				pnEdgePnts[nGp].dir = dir;
				pnEdgePnts[nGp].x = (int16)x;
				pnEdgePnts[nGp].y = (int16)y;
				pnEdgePnts[nGp].gx = pgx_calc[cPt];
				pnEdgePnts[nGp].gy = pgy_calc[cPt];
				nGp++;
			}
		}
	}
	time = clock() - time;

	//sorting the gradient points array
	qsort(pnEdgePnts, nGp, sizeof(pnEdgePnts[0]), sort_edgepnt);
	//printf("Angle values in [%d,%d]\n", pnEdgePnts[0].dir, pnEdgePnts[nGp-1].dir);
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
	int xB = -1, yB = -1;
	// main processing
	for (cGp = 0; cGp < nGp; ++cGp)
	{
		for (cPar = nBegPairIdx + cGp; cPar < nEndPairIdx + cGp; ++cPar)
		{
			D = pnEdgePnts[cPar].gx * pnEdgePnts[cGp].gy - pnEdgePnts[cPar].gy * pnEdgePnts[cGp].gx;
			if (D == 0)
				continue;

			Dt = pnEdgePnts[cPar].gx * (pnEdgePnts[cPar].y - pnEdgePnts[cGp].y) - pnEdgePnts[cPar].gy * (pnEdgePnts[cPar].x - pnEdgePnts[cGp].x);
			// t = (float)Dt / D;
			dx1 = pnEdgePnts[cGp].gx * Dt / D;
			dy1 = pnEdgePnts[cGp].gy * Dt / D;
			
			x = pnEdgePnts[cGp].x + dx1;
			y = pnEdgePnts[cGp].y + dy1;
			
			if (x < 1 || y < 1 || x > W - 1 || y > H - 1)
				continue;
			
			dx2 = pnEdgePnts[cPar].x - x;
			dy2 = pnEdgePnts[cPar].y - y;
			
			r1 = dx1*dx1 + dy1*dy1;
			r2 = dx2*dx2 + dy2*dy2;
			
			if (abs(r1 - r2) >= RADDELT || max(r1, r2) < MINRAD)
				continue;

			/*
			r1 = (r1 + r2) / 2;
			dx1 = pnEdgePnts[cGp].gx + pnEdgePnts[cPar].gx;
			dy1 = pnEdgePnts[cGp].gy + pnEdgePnts[cPar].gy;
			r2 = sqrt(dx1*dx1 + dx2*dx2);
			dx1 = (int)((float)dx1 / r2 * r1);
			dy1 = (int)((float)dy1 / r2 * r1);
			xB = max(min(W-10, x + dx1), 10);
			yB = max(min(H-10, y + dy1), 10);
			

			xB = max(min(W - 10, x + (int)r1), 10);
			yB = max(min(H - 10, y), 10);

			res = pucEdge[(yB-1)*W + xB-1] + pucEdge[(yB - 1) *W + xB] + pucEdge[(yB - 1)*W + xB+1] + 
				pucEdge[(yB + 0)*W + xB-1] + pucEdge[(yB + 0)*W + xB] + pucEdge[(yB + 0)*W + xB+1] + 
				pucEdge[(yB + 1)*W + xB-1] + pucEdge[(yB + 1)*W + xB] + pucEdge[(yB + 1)*W + xB+1];
			
			xB = max(min(W - 10, x - (int)r1), 10);
			yB = max(min(H - 10, y), 10);

			res = pucEdge[(yB - 1)*W + xB - 1] + pucEdge[(yB - 1) *W + xB] + pucEdge[(yB - 1)*W + xB + 1] +
				pucEdge[(yB + 0)*W + xB - 1] + pucEdge[(yB + 0)*W + xB] + pucEdge[(yB + 0)*W + xB + 1] +
				pucEdge[(yB + 1)*W + xB - 1] + pucEdge[(yB + 1)*W + xB] + pucEdge[(yB + 1)*W + xB + 1];
			pnAcc[y*W + x] += (res + 1);
			*/
			++pnAcc[y*W + x];
		}
	}
	time = clock() - time;
	printf("Paired Gradient voting: %f sec.\n", (float)time / CLOCKS_PER_SEC);

	maxRad = min(W, H) / 2;
	res = 0;
	if ((res = pgDetectCenter(&CI, &CP, pnAcc, H, W, flags)) != 0)
	{
		printf("ERROR: PG Accumulator processing: No center detected.");
		return -1;
	}
	
	vmax = rmax = 0;
	if (flags&BDL_PGRAD_CALCRADIUS)
	{
		/*
		if ((res = pgEstimateRadius(&CI, pnEdgePnts, nGp, H, W, flags)) != 0)
		{
			printf("[ERROR]: Radius estimation failed");
			return -1;
		}
		*/

		if ((res = pgEstimateRadii(&CP, &CI, pnEdgePnts, nGp, H, W, flags)) != 0)
		{
			printf("ERROR: Both radius estimation failed, code %d.", res);
			return -1;
		}
	}

	Result->IrisData->sCircle->r = CI.r;
	Result->IrisData->sCircle->xc = CI.xc;
	Result->IrisData->sCircle->yc = CI.yc;

	Result->PupilData->sCircle->r = CP.r;
	Result->PupilData->sCircle->xc = CP.xc;
	Result->PupilData->sCircle->yc = CP.yc;

	free(pnEdgePnts);
	free(pnAcc);
	free(pucEdge);
	free(canny_buf);
	return 0;
}

int CircBypass(int* radPath, uint8* priceMap)
{
	int i, n, sPnt, costUpd, c1, c2, c3, maxprice, optRoot, optRad, optCost;
	int *rootMap, *costMap, *parMap;
	
	const int H = NORMALIZED_IRIS_HEIGHT;
	const int W = CONTOUR_QUANTIZATION;
	rootMap = (int*)malloc(NORMALIZED_IRIS_HEIGHT * W * sizeof(int));
	costMap = (int*)malloc(NORMALIZED_IRIS_HEIGHT * W * sizeof(int));
	parMap = (int*)malloc(NORMALIZED_IRIS_HEIGHT * W * sizeof(int));


	maxprice = 2 * NORMALIZED_IRIS_HEIGHT * W;
	memset(costMap, maxprice, NORMALIZED_IRIS_HEIGHT * W * sizeof(int));
	memset(rootMap, -1, NORMALIZED_IRIS_HEIGHT * W * sizeof(int));
	memset(parMap, -1, NORMALIZED_IRIS_HEIGHT * W  * sizeof(int));

	for (i = 0; i < NORMALIZED_IRIS_HEIGHT; i++)
	{
		rootMap[i*W] = i;
		parMap[i*W] = i;
		costMap[i*W] = priceMap[i*W];
	}

	for (n = 1; n < W; n++)
	{
		c1 = costMap[n - 1] + (int)priceMap[n - 1];
		c2 = costMap[W + n - 1] + (int)priceMap[W + n - 1] + CDER;
		if (c1 < c2)
		{
			if (c1 >= 0)
			{
				costMap[n] = c1;
				rootMap[n] = rootMap[n - 1];
				parMap[n] = 0;
			}
			else
			{
				costMap[n] = c2;
				rootMap[n] = rootMap[n - 1 + W];
				parMap[n] = 1;
			}
		}
		else
		{
			if (c2 >= 0)
			{
				costMap[n] = c2;
				parMap[n] = 1;
				rootMap[n] = rootMap[W + n - 1];
			}
			else
			{
				costMap[n] = c1;
				rootMap[n] = rootMap[n - 1];
				parMap[n] = 0;
			}
		}

		for (i = 1; i < H - 1; i++)
		{
			c1 = costMap[i*W + n - 1] + (int)priceMap[i*W + n - 1];
			c2 = costMap[(i + 1)*W + n - 1] + (int)priceMap[(i + 1)*W + n - 1] + CDER;
			c3 = costMap[(i - 1)*W + n - 1] + (int)priceMap[(i - 1)*W + n - 1] + CDER;

			if (c1 < c2)
			{
				if (c1 < c3 && c1 >= 0)
				{
					costUpd = c1;
					sPnt = i;
				}
				else
				{
					if (c3 >= 0)
					{
						costUpd = c3;
						sPnt = i - 1;
					}
				}
			}
			else
			{
				if (c2 < c3)
				{
					if (c2 >= 0)
					{
						costUpd = c2;
						sPnt = i + 1;
					}
				}
				else
				{
					if (c3 >= 0)
					{
						costUpd = c3;
						sPnt = i - 1;
					}
				}
			}

			costMap[i*W + n] = costUpd;
			rootMap[i*W + n] = rootMap[sPnt * W + n - 1];
			parMap[i*W + n] = sPnt;
		}

		c1 = costMap[(H - 1) * W + n - 1] + (int)priceMap[(H - 1) * W + n - 1];
		c2 = costMap[(H - 2) * W + n - 1] + (int)priceMap[(H - 2) * W + n - 1] + CDER; //�������� ����� �������� ������ �� �������� �����������
		if (c1 < c2)
		{
			if (c1 >= 0)
			{
				costMap[(H - 1)*W + n] = c1;
				parMap[(H - 1)*W + n] = H - 1;
				rootMap[(H - 1)*W + n] = rootMap[(H - 1)*W + n - 1];
			}
			else
			{
				costMap[(H - 1)*W + n] = c2;
				parMap[(H - 1)*W + n] = H - 2;
				rootMap[(H - 1)*W + n] = rootMap[(H - 2)*W + n - 1];
			}
		}
		else
		{
			costMap[(H - 1)*W + n] = c2;
			parMap[(H - 1)*W + n] = H - 2;
			rootMap[(H - 1)*W + n] = rootMap[(H - 2)*W + n - 1];
		}
	}

	optCost = W * maxprice;
	for (sPnt = 0; sPnt < H; sPnt++)
	{
		if (rootMap[sPnt * W + W - 1] > sPnt + 2 || rootMap[sPnt * W + W - 1] < sPnt - 2)
			costMap[sPnt * W + W - 1] = maxprice;
		else
		{
			if (costMap[sPnt * W + W - 1] < optCost)
			{
				optCost = costMap[sPnt * W + W - 1];
				optRad = sPnt;
			}
		}
	}

	if (optCost > 0 && optCost <= CONTOUR_QUANTIZATION * maxprice)
	{
		optRoot = rootMap[optRad * W + W - 1];
		radPath[W] = (optRoot + optRad) / 2;
		radPath[W - 1] = optRad;
		for (i = W - 2; i >= 0; i--)
			radPath[i] = parMap[i + radPath[i + 1] * W];
	}

	free(rootMap);
	free(costMap);
	free(parMap);
	return optCost;
}


int IS_PupilCSP(SSegmentationResult *pResult, uint8 *edgeMap)
{
	int res, cost;
	int *radPath = (int*)malloc((CONTOUR_QUANTIZATION + 1) * sizeof(int));

	cost = CircBypass(radPath, edgeMap, pResult);
	if (cost < 0)
		return -1;

	pResult->PupilData->quality = 1.0 - (double)cost / (255 * CONTOUR_QUANTIZATION + CDER*CONTOUR_QUANTIZATION);

	int rMean = 0;
	int maxRad = (pResult->IrisData->sCircle->r + pResult->PupilData->sCircle->r) / 2;
	int minRad = 2 * pResult->PupilData->sCircle->r / 3;
	float factor = (float)(maxRad - minRad) / NORMALIZED_IRIS_HEIGHT;

	pResult->PupilData->sRCircle->xc = 0;
	pResult->PupilData->sRCircle->yc = 0;

	for (int i = 0; i < CONTOUR_QUANTIZATION; ++i)
	{
		float cRad = (float)(radPath[i])*factor + (float)minRad;
		rMean += cRad;

		edgeMap[i + radPath[i] * CONTOUR_QUANTIZATION] = 0;

		pResult->PupilData->sContour[i].x = pResult->PupilData->sCircle->xc + 
			(int)(cRad * cosf((float)i / CONTOUR_QUANTIZATION * 2 * PI));
		pResult->PupilData->sContour[i].y = pResult->PupilData->sCircle->yc +
			 (int)(cRad * sinf((float)i / CONTOUR_QUANTIZATION * 2 * PI));

		pResult->PupilData->sRCircle->xc += pResult->PupilData->sContour[i].x;
		pResult->PupilData->sRCircle->yc += pResult->PupilData->sContour[i].y;
	}
	
	SaveBmp8(pResult->name, "_PATH", CONTOUR_QUANTIZATION, NORMALIZED_IRIS_HEIGHT, edgeMap, 4);

	free(radPath);

	rMean /= CONTOUR_QUANTIZATION;
	pResult->PupilData->sRCircle->r = rMean;
	pResult->PupilData->sRCircle->xc /= CONTOUR_QUANTIZATION;
	pResult->PupilData->sRCircle->yc /= CONTOUR_QUANTIZATION;

	return 0;
}

int normalizeEdgeRing(SSegmentationResult *pResult, uint8 *normIris, uint8 *imgEdge, int Hc, int Wc, int xc, int yc)
{
	int maxRad = (pResult->IrisData->sCircle->r + pResult->PupilData->sCircle->r) / 2;
	int minRad = 2 * pResult->PupilData->sCircle->r / 3;
	double factor = (double)(maxRad - minRad) / NORMALIZED_IRIS_HEIGHT;
	int iCenter = Wc * yc + xc;

	for (int angle = 0; angle < CONTOUR_QUANTIZATION; ++angle)
	{
		for (int rad = 0; rad < NORMALIZED_IRIS_HEIGHT; ++rad)
		{
			int cRad = minRad + (int)(rad * factor);
			double x = (double)Wc/2 + (double)cRad * cos((double)(angle) / CONTOUR_QUANTIZATION * 2 * PI);
			double y = (double)Hc/2 + (double)cRad * sin((double)(angle) / CONTOUR_QUANTIZATION * 2 * PI);
				
			double interp = (1.0-frac(x))*(1.0-frac(y))*((int)imgEdge[(int)(y)*Wc + (int)(x)]) +
				(1.0 - frac(x))*frac(y)*((int)imgEdge[(int)ceil(y)*Wc + (int)(x)]) + 
				frac(x)*(1.0 - frac(y))*((int)imgEdge[(int)(y)*Wc + (int)ceil(x)]) + 
				frac(x)*frac(y)*((int)imgEdge[(int)ceil(y)*Wc + (int)ceil(x)]);
			normIris[rad * CONTOUR_QUANTIZATION + angle] = (uint8)interp;
		}
	}
	// SaveBmp8(pResult->name, "_NORMALIZED", CONTOUR_QUANTIZATION, NORMALIZED_IRIS_HEIGHT, normIris, 4);
	return 0;
}

int IS_RefinePupil(SSegmentationResult *pResult, const uint8* img, int H, int W, int flags)
{
	unsigned char *imgEdge, *imgCrop;
	int *imgGrad, *imgGradHist;
	void *canny_buf = NULL;
	int i, rad, angle, x, y, med, canny_sz, res;
	int gx, gy, ax, ay, max, perc;
	double t, tlow, thigh, sigma, val;
	int16 *pgx_calc, *pgy_calc;
	//FILE *fout = fopen("../data/res/edgeMapVals.txt", "w");

	//Cropping the pupil for better performance
	int xBeg = max(BORDER_VALUE, pResult->IrisData->sCircle->xc - pResult->IrisData->sCircle->r);
	int yBeg = max(BORDER_VALUE, pResult->IrisData->sCircle->yc - pResult->IrisData->sCircle->r);
	int xEnd = min(W - BORDER_VALUE, pResult->IrisData->sCircle->xc + pResult->IrisData->sCircle->r);
	int yEnd = min(W - BORDER_VALUE, pResult->IrisData->sCircle->yc + pResult->IrisData->sCircle->r);
	
	int Hc = yEnd - yBeg + 1;
	int Wc = xEnd - xBeg + 1;
	int cx = pResult->PupilData->sCircle->xc - xBeg;
	int cy = pResult->PupilData->sCircle->yc - yBeg;;

	if (Wc < W / 20 || Hc < H / 20)
		return -2;

	imgCrop = (uint8*)malloc(Hc * Wc * sizeof(uint8));
	imgGrad = (int*)malloc(Hc * Wc * sizeof(int));
	imgGradHist = (int*)malloc(Hc * Wc * sizeof(int));
	imgEdge = (unsigned char*)malloc(Hc * Wc * sizeof(unsigned char));

	memset(imgCrop, 0, Hc*Wc);
	memset(imgGrad, 0, Hc*Wc);
	canny_buf = (void*)malloc(20 * H * W * sizeof(unsigned char));

	for (y = yBeg; y <= yEnd; ++y)
	{
		for (x = xBeg; x <= xEnd; ++x)
		{
			imgCrop[(y - yBeg)* Wc + (x - xBeg)] = img[y * W + x];
		}
	}

	med = IPL_HIST_mdn_pixE(imgCrop, Hc*Wc);
	t = (double)med / 255;
	tlow = 0.8 * t;
	thigh = 1.33 * t;
	sigma = 2.f + 3 * t;
	canny_sz = 20 * H * W;//1049665; //�������
	if ((res = IPL_FILT_Canny(imgEdge, imgCrop, &pgx_calc, &pgy_calc, Wc, Hc, sigma, tlow, thigh,
		canny_buf, &canny_sz, "null")) != ERR_OK)					//original: 4.f .7f .8f
	{
		printf("Error %d: Canny, refinement stage fault\n", res);
		free(imgCrop);
		free(imgGrad);
		free(imgEdge);
		free(canny_buf);
		return -1;
	}

	int maxMagnitude = 0;
	for (y = 1; y < Hc; ++y)
	{
		for (x = 1; x < Wc - 1; ++x)
		{
			i = Wc * y + x;
			gx = (int)pgx_calc[i];
			gy = (int)pgy_calc[i];
			imgGrad[i] = (int)(sqrt(gx*gx + gy*gy));
			if (imgGrad[i] > maxMagnitude) maxMagnitude = imgGrad[i];
		}
	}
	memcpy(imgGradHist, imgGrad, Hc*Wc * sizeof(int));
	qsort(imgGradHist, Hc*Wc, sizeof(int), sort_int);
	perc = imgGradHist[Hc*Wc * 4 / 5];
	free(imgGradHist);
	for (y = 0; y < Hc; y++)
	{
		for (x = 0; x < Wc; x++)
		{
			i = Wc * y + x;
			ax = x - cx;
			ay = y - cy;
			gx = (int)pgx_calc[i];
			gy = (int)pgy_calc[i];
			val = (double)(ax * gx + ay * gy);
			val /= (int)(sqrt((float)(ax*ax + ay*ay)));
			val /= imgGrad[i];
			imgEdge[i] = (imgGrad > perc && val > 0.9) ? 255 - (uint8)((double)imgGrad[i] / maxMagnitude * 255) : 255;
		}
	}
	free(imgGrad);
	free(canny_buf);

	SaveBmp8(pResult->name, "_EDGE", Wc, Hc, imgEdge, 4);

	/*if (flags&BDL_PUPREF_SAVECROP)
		SaveBmp8("", "_cropcanny.bmp", Wc, Hc, imgEdge, 4);
	*/

	uint8 *normEdges = (uint8*)malloc(NORMALIZED_IRIS_HEIGHT * CONTOUR_QUANTIZATION * sizeof(uint8));
	memset(normEdges, 255, CONTOUR_QUANTIZATION * NORMALIZED_IRIS_HEIGHT * sizeof(uint8));

	if ((res = normalizeEdgeRing(pResult, normEdges, imgEdge, Hc, Wc, cx, cy)) != 0)
	{
		fprintf(stderr, "[ERROR] Iris edgemap normalization failed.");
		return -2;
	}

	free(imgEdge);
	free(imgCrop);

	if ((res = IS_PupilCSP(pResult, normEdges)) != 0)
	{
		fprintf(stderr, "Error %d: Circular shortest path obtaining failed: %d.\n", res);
		free(normEdges);
		return -3;
	}

	free(normEdges);
	printf("Contour quality: %f\n", pResult->PupilData->quality);

	return 0;
}