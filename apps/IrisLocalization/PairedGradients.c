#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "GlobalParams.h"
#include "PairedGradients.h"
#include "ImageProcessing.h"
#include "BorderDetection.h"

float frac(float val) {
	return val - floorf(val);
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


RESULT_CODE IS_Initialize(SSegmentationResult* pResult, char* name)
{
	if (pResult == NULL)
	{
		printf("\nError: Result struct was not initialized.");
		return ERROR_NULL_POINTER;
	}
	
	pResult->name = (char*)malloc(FILENAME_MAX * sizeof(char));
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
	pResult->PupilData->sContour = (sPnt*)malloc(CONTOUR_QUANTIZATION * sizeof(sPnt));
	pResult->PupilData->quality = -1.0;

	pResult->NormalizedIris = NULL; // Not defined what to do here yet

	return ERROR_OK;
}

RESULT_CODE IS_Deinitialize(SSegmentationResult* pResult)
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

	return ERROR_OK;
}


RESULT_CODE ImagePreprocessing(uint8* dstImg, uint8 *srcImg, int h, int w)
{
	// bad implementation
	IPL_FILT_HaussBlur3x3(dstImg, srcImg, w, h);
	
	uint8 *tmp = (uint8*)malloc(h*w * sizeof(uint8));
	if (NULL == tmp)
		return ERROR_NULL_POINTER;
	
	Dilate3x3Cross(tmp, dstImg, w, h);
	memcpy(dstImg, tmp, h*w * sizeof(uint8));
	free(tmp);
	return ERROR_OK;
}


RESULT_CODE IrisSegmentation(SSegmentationResult* Result, //Segmentation result
	uint8 *imgInput, int H, int W,
	int angle, int flags)
{
	if (Result == NULL)
	{
		printf("\n[ERROR]: IrisSegmentation: Segmentation result was not initialized.");
		return -1;
	}
	int res;
	uint8* dstImg;
	if (flags | BDL_PGRAD_PREPROCESSING)
	{	
		dstImg = (uint8*)malloc(H*W*sizeof(uint8));
		// pre-processing stage
		if ((res = ImagePreprocessing(dstImg, imgInput, H, W)) < 0)
		{
			free(dstImg);
			fprintf(stderr, "\n[ERROR] IrisSegmentation(): Preprocessing fail.");
			return res;
		}
	}
	else
		dstImg = imgInput;
	
	//Apply the gradient pairs method to the image to find the first boundary
	if ((res = IS_ApproxBothBorders(Result, dstImg, H, W, angle, flags)) != 0)
	{	
		if (flags | BDL_PGRAD_PREPROCESSING)
			free(dstImg);
		fprintf(stderr, "[ERROR] IrisSegmentation(): Gradient pair method fail.\n");
		return res;
	}
	
	//Apply the circular shortest path refinement method to the found pupil circle approximation
	if ((res = IS_RefinePupil(Result, dstImg, H, W, flags)) != 0)
	{
		if (flags | BDL_PGRAD_PREPROCESSING)
			free(dstImg);
		fprintf(stderr, "Error: Refinement fail\n");
		return res;
	}
	
	if (flags | BDL_PGRAD_PREPROCESSING)
		free(dstImg);
		
	return 0;
}

/* Returns binarization threshold for grayscale image with size pixels. */
uint8 otsuThreshold(const unsigned char *image, const int size)
{
	uint8 min = image[0], max = image[0];
	uint8 tmp;
	int temp, temp1;
	int *hist;
	int histSize;

	int alpha, beta, threshold = 0;
	float sigma, maxSigma = -1;
	float w1, a;

	/**** Compute the histogram ****/
	/* Search for the largest and the smallest semitones */
	for (int i = 1; i < size; ++i)
	{
		tmp = image[i];
		if (tmp < min)   min = tmp;
		if (tmp > max)   max = tmp;
	}

	histSize = (int)(max - min) + 1;
	if ((hist = (int*) malloc(sizeof(int) * histSize)) == NULL) return ERROR_NULL_POINTER;

	memset(hist, 0, histSize * sizeof(int));

	/* Calculate semitones */
	for (int i = 0; i < size; ++i)
		++hist[(int)(image[i] - min)];

	/**** Histogram is ready ****/

	temp = temp1 = 0;
	alpha = beta = 0;
	/* Calculating Expectations */
	for (int i = 0; i <= (int)(max - min); ++i)
	{
		temp += i*hist[i];
		temp1 += hist[i];
	}

	/* Main cycle for threshold calculation
	Cycle through all the semi-tones to search for the one with the smallest intra-clss variation */
	for (int i = 0; i<(int)(max - min); ++i)
	{
		alpha += i*hist[i];
		beta += hist[i];

		w1 = (float)beta / temp1;
		a = (float)alpha / beta - (float)(temp - alpha) / (temp1 - beta);
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
RESULT_CODE pgDetectCenter(SCircleData *CI, SCircleData *CP, int* pnAcc, int H, int W, int flags)
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
				pnAccOut[y*W + x] = (uint8)((float)pnAccBl[y*W + x] / vmax * 255);
		/*
		SaveBmp8(name, "_acc1.bmp", W, H, pnAccOut, 4);
		*/
		free(pnAccOut);
	}
	
	free(pnAccBl);

	if (vmax == 0)
		return ERROR_NO_ANSWER;

	return ERROR_OK;
}


RESULT_CODE pgEstimateRadii(SCircleData *CP, SCircleData *CI, SEdgePnt* pnEdgePnts, int nGp, int H, int W, int flags)
{
	int dx1, dy1, gx, gy, r, rMax, vrmax, vmax, cGp;
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
		radHistSmooth[r] = 0;

	hRmax = 4 * vmax / 3;
	int lRmax = 3 * vmax / 4;
	for (r = lRmax; r < hRmax; ++r)
		radHistSmooth[r] = 0;

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
		return ERROR_NO_ANSWER;

	return ERROR_OK;
}


RESULT_CODE GetEdgePointsFromSobel(uint8* pucEdge, int16** pGx, int16** pGy,
	const uint8* img, int W, int H, int nPercHigh, int nPercLow, int buffer_size)
{
	return ERROR_NO_ANSWER;
}



RESULT_CODE GetEdgePointsFromCanny(uint8* pucEdge, int16** pGx, int16** pGy,
	const uint8* img, int W, int H, float tHigh, float tLow, float sigma, int buffer_size)
{
	int canny_sz, res;
	void* canny_buf = NULL;
	canny_buf = (void*)malloc(4 * buffer_size);

	//Canny edge detection and intensity gradients calculating	
	canny_sz = 4 * buffer_size; //5222400
	if ((res = IPL_FILT_Canny(pucEdge, img, pGx, pGy, W, H,
		sigma, //5.f,
		tLow,
		tHigh,
		canny_buf, &canny_sz, "name")) != ERROR_OK)					//original: 4.f .7f .8f
	{
		printf("[ERROR]: Canny failed %d.\n", res);
		return -1;
	}
	return ERROR_OK;
}

void GetCannyParamsFromOtsu(float* tHigh, float* tLow, float* sigma, const uint8* img, int H, int W)
{
	//median search for canny tresholding
	uint8 uMax;
	uint8* imgsort = (uint8*)malloc(H*W * sizeof(imgsort[0]));
	memcpy(imgsort, img, sizeof(img[0])*H*W);
	//TODO: Estimate image median fast
	qsort(imgsort, H*W, sizeof(uint8), sort_img8U);

	uMax = imgsort[H*W - 1];
	free(imgsort);

	uint8 ucOtsu = otsuThreshold(img, H*W);
	/*
	uint8 *imgOtsu = (uint8*)malloc(H*W * sizeof(uint8));
	for (int i = 0; i < H * W; ++i)
		imgOtsu[i] = img[i] > ucOtsu ? 0 : 255;

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
	free(imgOtsu);

	for (int x = 0; x < W; ++x)
	{
		xProj[x] /= H;
		if (xProj != 255)
			break;
		xBeg = x;
	}
	for (int y = 0; y < H; ++y)
		yProj[y] /= W;
	*/

	// printf("[PG] t_otsu = %f\n", (float)ucOtsu / 255);
	float tUpper = 1.2f * (float)ucOtsu / (H * W);
	*tHigh = MIA_min(.7f, tUpper);
	*tLow = *tHigh / 2;
	//tlow = .55f;
	//thigh = .7f;
	*sigma = 4.f;
}

#define LOG_EM_PG	0
RESULT_CODE DetectCenterFromEdgeMapPairedGradients(int* cenX, int* cenY, SEdgePnt** ppnEdgePnts, int* nGp,
	const uint8* img, int H, int W, uint8* pucEdge, int16* pGx,int16* pGy, int dPhi, int flags)
{
	int res = ERROR_OK;
	int *pnAcc = NULL;

	int sz = 0;
	sz += H*W * sizeof(*ppnEdgePnts[0]) +  // edge array
		H*W * sizeof(pnAcc[0]);     // accum
									//H*W * sizeof(pnCircleCnd[0]); // circle params 
	
	SEdgePnt* pnEdgePnts;
	short toFreeEdgePntArrray = 0;
	if (NULL == ppnEdgePnts)
	{
		pnEdgePnts = (SEdgePnt*)malloc(H*W * sizeof(pnEdgePnts[0]));
#if LOG_EM_PG
		printf("[LOG] DetectCenterFromEdgeMapPairedGradients(): Allocated new memory for Edge Points list.\n");
#endif
		toFreeEdgePntArrray = 1;
	}
	else
	{
		pnEdgePnts = *ppnEdgePnts;
	}

	pnAcc = (int*)malloc(H*W * sizeof(pnAcc[0]));
#if LOG_EM_PG
	printf("[LOG] DetectCenterFromEdgeMapPairedGradients(): Allocated new memory for Accumulator.\n");
#endif
#if TIMERS_ENABLED
	time_t time = clock();
#endif
	// Collect gradient point to array
	*nGp = 0;
	int x, y, cPt, dir;
	for (y = 10; y < H - 10; ++y)
	{
		for (x = 10; x < W - 10; ++x)
		{
			if (pucEdge[cPt = y*W + x])
			{
				dir = (int)(.5 + 180.*atan2f((float)pGy[cPt], (float)pGx[cPt]) / fPI);
				pnEdgePnts[*nGp].dir = dir;
				pnEdgePnts[*nGp].x = (int16)x;
				pnEdgePnts[*nGp].y = (int16)y;
				pnEdgePnts[*nGp].gx = pGx[cPt];
				pnEdgePnts[*nGp].gy = pGy[cPt];
				++(*nGp);
			}
		}
	}
#if TIMERS_ENABLED
	time = clock() - time;
#endif
#if LOG_EM_PG
	printf("[LOG] DetectCenterFromEdgeMapPairedGradients(): Filled EdgePoints list: %d entries.\n", *nGp);
#endif

	//sorting the gradient points array
	qsort(pnEdgePnts, *nGp, sizeof(pnEdgePnts[0]), sort_edgepnt);
	//duplicate the array
	memcpy(pnEdgePnts + *nGp, pnEdgePnts, sizeof(pnEdgePnts[0]) * *nGp);
	//clean the accumulator
	for (int cGp = 0; cGp<*nGp; cGp++)
		pnEdgePnts[cGp + *nGp].dir = pnEdgePnts[cGp].dir + 360;
	// clean separated accumulators for center and radius
	memset(pnAcc, 0, W*H * sizeof(pnAcc[0]));

	// prime counters
	int nBegPairIdx = 0, nEndPairIdx = 0;
	int tint = pnEdgePnts[0].dir + dPhi - ANGDELT;
	for (nBegPairIdx = 0; ((nBegPairIdx<*nGp) &&
		(pnEdgePnts[nBegPairIdx].dir<tint)); nBegPairIdx++);

	tint = pnEdgePnts[0].dir + dPhi + ANGDELT;
	for (nEndPairIdx = nBegPairIdx; ((nEndPairIdx<*nGp) &&
		(pnEdgePnts[nEndPairIdx].dir<tint)); nEndPairIdx++);

	int nC = 0;
#if TIMERS_ENABLED
	time = clock();
#endif
	int xB = -1, yB = -1;
	int xmax = -1, ymax = -1, rmax = -1;
	int dx1, dx2, dy1, dy2, r1, r2;
	// main processing
	for (int cGp = 0; cGp < *nGp; ++cGp)
	{
		for (int cPar = nBegPairIdx + cGp; cPar < nEndPairIdx + cGp; ++cPar)
		{
			int D = pnEdgePnts[cPar].gx * pnEdgePnts[cGp].gy - pnEdgePnts[cPar].gy * pnEdgePnts[cGp].gx;
			if (D == 0)
				continue;

			int Dt = pnEdgePnts[cPar].gx * (pnEdgePnts[cPar].y - pnEdgePnts[cGp].y) - pnEdgePnts[cPar].gy * (pnEdgePnts[cPar].x - pnEdgePnts[cGp].x);
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
#if LOG_EM_PG
	printf("[LOG] DetectCenterFromEdgeMapPairedGradients(): Paired Gradient voting: finished.\n");
#endif
#if TIMERS_ENABLED
	time = clock() - time;
	printf("Paired Gradient voting: %f sec.\n", (float)time / CLOCKS_PER_SEC);
#endif
	SCircleData CI, CP;
	res = pgDetectCenter(&CI, &CP, pnAcc, H, W, flags);
	if (ERROR_OK != res)
	{
		printf("[ERROR] DetectCenterFromEdgeMapPairedGradients(): PG Accumulator processing: No center detected.\n");
		return res;
	}
	*cenX = CI.xc;
	*cenY = CI.yc;

	if (toFreeEdgePntArrray)
		free(pnEdgePnts);

	free(pnAcc);
	return ERROR_OK;

}

RESULT_CODE DetectCenterPairedGradients(int* cenX, int* cenY, SEdgePnt** ppnEdgePnts, int* nGp,
	const uint8* img, int H, int W, float tHigh, float tLow, float sigma, int dPhi, int flags)
{
	int res = ERROR_OK;
	int *pnAcc = NULL;

	int sz = 0;
	sz += H*W * sizeof(*ppnEdgePnts[0]) +  // edge array
		H*W * sizeof(pnAcc[0]);     // accum
									//H*W * sizeof(pnCircleCnd[0]); // circle params 
	int16 *pGx, *pGy;
	uint8* pucEdge = (uint8*)malloc(H*W * sizeof(pucEdge[0]));
	if (flags & BDL_PGRAD_USE_CANNY)
	{
		float _tHigh, _tLow, _sigma;
		if (flags & BDL_PGRAD_CANNY_AUTO)
		{
			GetCannyParamsFromOtsu(&_tHigh, &_tLow, &_sigma, img, H, W);
			// printf("[LOG] DetectCenterPairedGradients(): used automatic Canny params set up.");
		}
		else
		{
			_tHigh = tHigh;
			_tLow = tLow;
			_sigma = sigma;
			// printf("[LOG] DetectCenterPairedGradients(): used Canny params from function input.");
		}

		// printf("LOG: Got Canny params from otsu: %f, %f, %f", tHigh, tLow, sigma);
		res = GetEdgePointsFromCanny(pucEdge, &pGx, &pGy, img, W, H, _tHigh, _tLow, _sigma, 4 * sz);
		// printf("LOG: Got Edge points from Canny");

		if (ERROR_OK != res)
		{
			// printf("[Error]: Canny edge detection failed with result code %d.\n", res);
			free(pucEdge);
			return res;
		}
	}
	else if (flags & BDL_PGRAD_USE_SOBEL)
	{
		// not implemented yet
		return ERROR_WRONG_INPUT;
	}
	else if (flags & BDL_PGRAD_USE_IMAGE)
	{
		pucEdge = (uint8*)img;
		pGx = (int16*)(pucEdge + H*W * sizeof(int16*));
		pGy = (int16*)(pGx + H*W * sizeof(int16*));
	}
	else
	{
		// not implemented yet
		return ERROR_WRONG_INPUT;
	}
	res = DetectCenterFromEdgeMapPairedGradients(cenX, cenY, ppnEdgePnts, nGp, img, H, W, pucEdge, pGx, pGy, dPhi, flags);
	free(pucEdge);
	return res;
}

RESULT_CODE IS_ApproxBothBorders(SSegmentationResult *Result, const uint8 *img, int H, int W, int dAngle, int flags)
{
#if TIMERS_ENABLED
	clock_t time = clock();
#endif
	SCircleData CI, CP;
	int res, cenX, cenY;
	
	// Detect center with paired gradients
	SEdgePnt* pnEdgePnts = (SEdgePnt*)malloc(H*W*sizeof(SEdgePnt));
	int nGp = 0;
	res = DetectCenterPairedGradients(&cenX, &cenY, &pnEdgePnts, &nGp, img, H, W, .4f, .7f, 4.0f, dAngle, flags);
	if (ERROR_OK != res)
	{
		printf("[Error] Center detection with paired gradients failed with error code %d.", res);
		return res;
	}
#if TIMERS_ENABLED
	time = clock() - time;
	printf("Paired Gradient voting: %f sec.\n", (float)time / CLOCKS_PER_SEC);
#endif
	
	int vmax = 0, rmax = 0;
	res = pgEstimateRadii(&CP, &CI, pnEdgePnts, nGp, H, W, flags);
	if (ERROR_OK != res)
	{
		printf("ERROR: Both radius estimation failed, code %d.", res);
		return res;
	}
	free(pnEdgePnts);
	Result->IrisData->sCircle->r = CI.r;
	Result->IrisData->sCircle->xc = CI.xc;
	Result->IrisData->sCircle->yc = CI.yc;

	Result->PupilData->sCircle->r = CP.r;
	Result->PupilData->sCircle->xc = CP.xc;
	Result->PupilData->sCircle->yc = CP.yc;
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

	// initialization
	maxprice = 2 * NORMALIZED_IRIS_HEIGHT * W;
	for (i = 0; i < NORMALIZED_IRIS_HEIGHT * W; ++i)
	{
		costMap[i] = maxprice;
		rootMap[i] = -1;
		parMap[i] = -1;
	}

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
		c2 = costMap[(H - 2) * W + n - 1] + (int)priceMap[(H - 2) * W + n - 1] + CDER; //Изменить потом величину штрафа за величину производной
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


RESULT_CODE IS_PupilCSP(SSegmentationResult *pResult, uint8 *edgeMap)
{
	int cost;
	int *radPath = (int*)malloc((CONTOUR_QUANTIZATION + 1) * sizeof(int));

	cost = CircBypass(radPath, edgeMap);//, pResult);
	if (cost < 0)
		return ERROR_NO_ANSWER;

	pResult->PupilData->quality = 1.0f - (float)cost / (255 * CONTOUR_QUANTIZATION + CDER*CONTOUR_QUANTIZATION);

	int rMean = 0;
	int maxRad = (pResult->IrisData->sCircle->r + pResult->PupilData->sCircle->r) / 2;
	int minRad = 2 * pResult->PupilData->sCircle->r / 3;
	float factor = (float)(maxRad - minRad) / NORMALIZED_IRIS_HEIGHT;

	pResult->PupilData->sRCircle->xc = 0;
	pResult->PupilData->sRCircle->yc = 0;

	for (int i = 0; i < CONTOUR_QUANTIZATION; ++i)
	{
		float cRad = (float)(radPath[i])*factor + (float)minRad;
		rMean += (int)(cRad+0.5f);

		edgeMap[i + radPath[i] * CONTOUR_QUANTIZATION] = 0;

		pResult->PupilData->sContour[i].x = pResult->PupilData->sCircle->xc + 
			(int)(cRad * cosf((float)i / CONTOUR_QUANTIZATION * 2 * PI));
		pResult->PupilData->sContour[i].y = pResult->PupilData->sCircle->yc +
			 (int)(cRad * sinf((float)i / CONTOUR_QUANTIZATION * 2 * PI));

		pResult->PupilData->sRCircle->xc += pResult->PupilData->sContour[i].x;
		pResult->PupilData->sRCircle->yc += pResult->PupilData->sContour[i].y;
	}

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
			float x = (float)Wc/2 + (float)cRad * cosf((float)(angle) / CONTOUR_QUANTIZATION * 2 * fPI);
			float y = (float)Hc/2 + (float)cRad * sinf((float)(angle) / CONTOUR_QUANTIZATION * 2 * fPI);
				
			float interp = (1.0f-frac(x))*(1.0f-frac(y))*((int)imgEdge[(int)(y)*Wc + (int)(x)]) +
				(1.0f - frac(x))*frac(y)*((int)imgEdge[(int)ceilf(y)*Wc + (int)(x)]) + 
				frac(x)*(1.0f - frac(y))*((int)imgEdge[(int)(y)*Wc + (int)ceilf(x)]) + 
				frac(x)*frac(y)*((int)imgEdge[(int)ceilf(y)*Wc + (int)ceilf(x)]);
			normIris[rad * CONTOUR_QUANTIZATION + angle] = (uint8)interp;
		}
	}
	// SaveBmp8(pResult->name, "_NORMALIZED", CONTOUR_QUANTIZATION, NORMALIZED_IRIS_HEIGHT, normIris, 4);
	return 0;
}

RESULT_CODE IS_RefinePupil(SSegmentationResult *pResult, const uint8* img, int H, int W, int flags)
{
	unsigned char *imgEdge, *imgCrop;
	int *imgGrad, *imgGradHist;
	void *canny_buf = NULL;
	int i, x, y, med, canny_sz, res;
	int gx, gy, ax, ay,  perc;
	double t, tlow, thigh, sigma;
	float val;
	int16 *pgx_calc, *pgy_calc;

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
	canny_sz = 20 * H * W;//1049665; // Crutch
	if ((res = IPL_FILT_Canny(imgEdge, imgCrop, &pgx_calc, &pgy_calc, Wc, Hc, sigma, tlow, thigh,
		canny_buf, &canny_sz, "null")) != ERROR_OK)					//original: 4.f .7f .8f
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
			val = (float)(ax * gx + ay * gy);
			val /= (int)(sqrt((float)(ax*ax + ay*ay)));
			val /= imgGrad[i];
			imgEdge[i] = (imgGrad[i] > perc && val > 0.9f) ? 255 - (uint8)((float)imgGrad[i] / maxMagnitude * 255) : 255;
		}
	}
	free(imgGrad);
	free(canny_buf);

	// SaveBmp8(pResult->name, "_EDGE", Wc, Hc, imgEdge, 4);

	/*if (flags&BDL_PUPREF_SAVECROP)
		SaveBmp8("", "_cropcanny.bmp", Wc, Hc, imgEdge, 4);
	*/

	uint8 *normEdges = (uint8*)malloc(NORMALIZED_IRIS_HEIGHT * CONTOUR_QUANTIZATION * sizeof(uint8));
	memset(normEdges, 255, CONTOUR_QUANTIZATION * NORMALIZED_IRIS_HEIGHT * sizeof(uint8));

	if ((res = normalizeEdgeRing(pResult, normEdges, imgEdge, Hc, Wc, cx, cy)) != 0)
	{
		fprintf(stderr, "[ERROR] Iris edgemap normalization failed: error %d.\n", res);
		return -2;
	}

	free(imgEdge);
	free(imgCrop);

	if ((res = IS_PupilCSP(pResult, normEdges)) != 0)
	{
		fprintf(stderr, "[Error]: Circular shortest path obtaining failed: %d.\n", res);
		free(normEdges);
		return -3;
	}

	free(normEdges);
	printf("Contour quality: %f\n", pResult->PupilData->quality);

	return 0;
}