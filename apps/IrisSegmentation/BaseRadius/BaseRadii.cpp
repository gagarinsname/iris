/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  11 November 2017                                       */
/*      Purpose:  Iris library API                                       */
/*      Authors:  Iurii Efimov                                           */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <iostream>

#include "GlobalParams.h"
#include "BaseRadii.h"

#define INFINITY FLT_MAX / 10 // max_float

float metricGradWindowDTW(int* pDataL, int posL, int* pDataR, int posR, int hwSize)
{
	float gradRes = .0f, gLeft, gRight;
	for (int i = -hwSize + 2; i <= hwSize - 2; ++i)
	{
		// gradient
		gLeft = -1 * pDataL[posL + i - 2] - 2 * pDataL[posL + i - 1] + 2 * pDataL[posL + i + 1] + 1 * pDataL[posL + i + 2];
		gRight = -1 * pDataR[posR + i - 2] - 2 * pDataR[posR + i - 1] + 2 * pDataR[posR + i + 1] + 1 * pDataR[posR + i + 2];
		gradRes += (float)((gLeft - gRight)*(gLeft - gRight));
	}
	return 1.0;
}

float metricWindowDTW(int* pDataL, int posL, int* pDataR, int posR, int hwSize)
{
	float res = 0.0f;
	for (int i = -hwSize; i <= hwSize; ++i)
	{
		// euclidean
		// fprintf(stdout, "%d %d %d %d\n", posR, posR + i, posL, posL + i);
		res += (float)((pDataL[posL + i] - pDataR[posR + i])*(pDataL[posL + i] - pDataR[posR + i]));
	}
	res = sqrtf(res) / (2 * hwSize + 1);
	// fprintf(stdout, "metric res: %f\n", res);
	return res;
}

void downsampleSeries(int* dst, int dstSize, int* src, int srcSize, int factor)
{
	int i;
	if (1 == factor)
	{
		memcpy(dst, src, srcSize * sizeof(int));
		return;
	}
	// memset(dst, 0, dstSize * sizeof(int));
	int* temp = (int*)malloc(srcSize * sizeof(int));
	int hw = factor / 2;
	IPL_HIST_Blur(temp, (const int*)src, srcSize, hw);

	for (i = 0; i < dstSize; ++i)
	{
		dst[i] = temp[factor * i] / (2 * hw);
	}
	free(temp);
}

void upsamplePath(sPnt* path, int* sizeP, int scaleFactor)
{
	sPnt* tmp = (sPnt*)malloc(*sizeP * sizeof(sPnt));
	for (int i = 0; i < *sizeP; ++i) {
		tmp[i].x = path[i].x;
		tmp[i].y = path[i].y;
	}

	for (int i = 0; i < *sizeP; ++i)
	{
		int k = scaleFactor * i;
		for (int j = 0; j < scaleFactor; ++j) {
			path[k + j].x = scaleFactor * tmp[i].x;
			path[k + j].y = scaleFactor * tmp[i].y + j;
		}
	}
	*sizeP = *sizeP * scaleFactor;
	free(tmp);
}

RESULT_CODE sequenceGradient(int *pGrData, int *pData, int size)
{
	const int kSize = 3;
	int kernel[] = { -1,-4,-6,0,6,4,1 }; // convolving 1x5 gaussian kernel with [-1,0,1]
	memset(pGrData, 0, size * sizeof(int));
	if (size < 2 * kSize + 1)
	{
		return ERROR_WRONG_INPUT;
	}
	// padding: ignoring edge regions
	// main part
	int* cData = &(pData[kSize]);
	for (int i = kSize; i < size - kSize; ++i)
	{
		pGrData[i] = 0;
		for (int j = -kSize; j < kSize; ++j)
			pGrData[i] += kernel[kSize+j] * cData[-j];
		++cData;
	}
	return ERROR_OK;
}




RESULT_CODE calcFullDynamicTimeWarping(sPnt* pathDTW, int* sizePath, float* matDTW,
	int* pDataL, int* pGradL, int sizeL, int* pDataR, int* pGradR, int sizeR,
	int compRegion, float regWeight)
{
	// Previous state matrix for minimal cost extraction
	sPnt* matrixPrevState = (sPnt*)malloc((sizeL)* (sizeR) * sizeof(sPnt));
	int res = ERROR_OK;
	
	matDTW[0] = 0.0f;
	for (int i = 0; i < sizeL; ++i)
	{
		for (int j = 0; j < sizeR; ++j)
		{
			matDTW[i * sizeR + j] = INFINITY;
			matrixPrevState[i * sizeR + j].x = (i < j) ? j - 1 : j;
			matrixPrevState[i * sizeR + j].y = (i < j) ? i : i - 1;
		}
	}
	// correct distance padding
	for (int k = 0; k < sizeL; ++k)
	{
		int dValue = pDataR[0] - pDataL[k], dGrad = (int)(pGradR[0] > 0 ? 1 : -1) - (int)(pGradL[k] > 0 ? 1 : -1);;
		matDTW[k*sizeR] = (float)(MIA_abs(dValue));// + MIA_abs(dGrad);
		matDTW[k] = (float)(MIA_abs(pDataR[k] - pDataL[0])); // +MIA_abs(pGradR[k] - pGradL[0]);
	}

	if (1 > compRegion)
	{
		for (int i = 1; i < sizeL; ++i)
		{
			int jBeg = MIA_max(1, i - 20), jEnd = MIA_min(sizeR, i + 20);
			for (int j = jBeg; j < jEnd; ++j)
			{
				float* ptrC = &matDTW[i* sizeR + j];
				float* ptrT = &matDTW[(i - 1) * sizeR + j];
				float* ptrL = &matDTW[i * sizeR + j - 1];
				float* ptrTL = &matDTW[(i - 1) * sizeR + j - 1];
				sPnt* ptrPrev = &matrixPrevState[i * sizeR + j];

				// *ptrC = MIA_abs(pDataR[j] - pDataL[i]) + MIA_abs(pGradR[j] - pGradL[i]);
				int dValue = pDataR[j] - pDataL[i], dGrad = (int)(pGradR[j] > 0 ? 1 : -1) - (int)(pGradL[i] > 0 ? 1 : -1);
				int dValue_0 = pDataR[j-1] - pDataL[i-1], dGrad_0 = pGradR[j-1] - pGradL[i-1];
				*ptrC = sqrtf((float)(dValue*dValue + dValue_0*dValue_0)) + regWeight * (float)((i - j) * (i - j)) + (float)(MIA_abs(dGrad));
				float tmp = *ptrC;

				if (*ptrTL <= MIA_min(*ptrL, *ptrT))
				{
					
					*ptrC = tmp + *ptrTL;
					ptrPrev->x = j - 1;
					ptrPrev->y = i - 1;
				}
				else if (*ptrL <= MIA_min(*ptrT, *ptrTL))
				{
					*ptrC = tmp + *ptrL;
					ptrPrev->x = j - 1;
					ptrPrev->y = i;
				}
				else
				{
					*ptrC = tmp + *ptrT;
					ptrPrev->x = j;
					ptrPrev->y = i - 1;
				}
				// ++ptrC; ++ptrL; ++ptrT; ++ptrTL;

			}
		}
		// search for the minimal total cost
		float optPathValue = matDTW[(sizeL - 1) * sizeR + sizeR - 1];
		sPnt optPathEnd = { sizeR - 1, sizeL - 1 };
		// obtain minimal total cost path
		for (int i = 0; i < *sizePath; ++i)
		{
			int xEnd = optPathEnd.x, yEnd = optPathEnd.y;
			// std::cout << "i: " << i << " xEnd: " << xEnd << " yEnd: " << yEnd << std::endl;
			if (xEnd < 0 || yEnd < 0 || xEnd >= sizeR || yEnd >= sizeL)
			{
				*sizePath = i;
				break;
			}
			pathDTW[i].x = optPathEnd.x;
			pathDTW[i].y = optPathEnd.y;
			optPathEnd.x = matrixPrevState[yEnd * sizeR + xEnd].x;
			optPathEnd.y = matrixPrevState[yEnd * sizeR + xEnd].y;
		}
		for (int i = 0; i < *sizePath / 2; ++i)
		{
			sPnt t;
			t.x = pathDTW[i].x;
			t.y = pathDTW[i].y;
			pathDTW[i].x = pathDTW[*sizePath - i - 1].x;
			pathDTW[i].y = pathDTW[*sizePath - i - 1].y;
			pathDTW[*sizePath - i - 1].x = t.x;
			pathDTW[*sizePath - i - 1].y = t.y;
		}
	}
	else
	{
		for (int k = 0; k < *sizePath; ++k)
		{
			int i = pathDTW[k].y;
			if (i == 0)
				continue;
			int jBeg = MIA_max(1, pathDTW[k].x - compRegion);
			int jEnd = MIA_min(sizeR, pathDTW[k].x + compRegion);
			
			float* ptrC = matDTW + (i)* sizeR + jBeg;
			if (*ptrC < INFINITY) continue;
			std::cout << "LOG: i:" << i << " jRange:  jBeg: " << jBeg << ", jEnd: " << jEnd << std::endl;

			float* ptrT = matDTW + (i - 1) * sizeR + jBeg;
			float* ptrL = matDTW + i * sizeR + jBeg - 1;
			float* ptrTL = matDTW + (i - 1) * sizeR + jBeg - 1;
			sPnt* ptrPrev = matrixPrevState + i * sizeR + jBeg;
			
			for (int j = jBeg; j < jEnd; ++j)
			{
				// *ptrC = MIA_abs(pDataR[j] - pDataL[i]) + MIA_abs(pGradR[j] - pGradL[i]);
				int dValue = pDataR[j] - pDataL[i], dGrad = pGradR[j] - pGradL[i];
				*ptrC = INTSQRT(dValue*dValue + dGrad*dGrad);

				if (*ptrTL <= MIA_min(*ptrL, *ptrT))
				{
					*ptrC += *ptrTL;
					ptrPrev->x = j - 1;
					ptrPrev->y = i - 1;
				}
				else if (*ptrL <= MIA_min(*ptrT, *ptrTL))
				{
					*ptrC += *ptrL;
					ptrPrev->x = j - 1;
					ptrPrev->y = i;
				}
				else
				{
					*ptrC += *ptrT;
					ptrPrev->x = j;
					ptrPrev->y = i - 1;
				}
				++ptrPrev; ++ptrC; ++ptrL; ++ptrT; ++ptrTL;

			}
		}
		// search for the minimal total cost
		float optPathValue = matDTW[sizeL* sizeR - 1];
		sPnt optPathEnd = { sizeR - 1, sizeL - 1 };
		// float* ptrLastCol = matDTW + sizeR * sizeL - 1;
		// float* ptrLastRow = ptrLastCol;
		// obtain minimal total cost path
		for (int i = 0; i < *sizePath; ++i)
		{
			int xEnd = optPathEnd.x, yEnd = optPathEnd.y;
			if (xEnd < 0 || yEnd < 0)
			{
				*sizePath = i;
				std::cout << "Warning: path is out of range: x=" << xEnd << " y=" << yEnd << std::endl;
				break;
			}
			
			pathDTW[i].x = optPathEnd.x;
			pathDTW[i].y = optPathEnd.y;

			optPathEnd.x = matrixPrevState[yEnd * sizeR + xEnd].x;
			optPathEnd.y = matrixPrevState[yEnd * sizeR + xEnd].y;
		}
		for (int i = 0; i < *sizePath / 2; ++i)
		{
			sPnt t;
			t.x = pathDTW[i].x;
			t.y = pathDTW[i].y;
			pathDTW[i].x = pathDTW[*sizePath - i - 1].x;
			pathDTW[i].y = pathDTW[*sizePath - i - 1].y;
			pathDTW[*sizePath - i - 1].x = t.x;
			pathDTW[*sizePath - i - 1].y = t.y;
		}
	}
	free(matrixPrevState);
	return res;
}

RESULT_CODE calcFullDynamicTimeWarpingWrapper(sPnt* pathDTW, int* sizePath, float* matDTW, int* pDataL, int sizeL, int* pDataR, int sizeR,
	int nCompare, float regWeight)
{
	if (sizeR <= 0 || sizeL <= 0 || NULL == pDataR || NULL == pDataL)
		return ERROR_WRONG_INPUT;

	int* pNGradR = (int*)malloc(sizeR * sizeof(int));
	sequenceGradient(pNGradR, pDataR, sizeR);

	int* pNGradL = (int*)malloc(sizeL * sizeof(int));
	sequenceGradient(pNGradL, pDataL, sizeL);

	int res = calcFullDynamicTimeWarping(pathDTW, sizePath, matDTW, pDataL, pNGradL, sizeL, pDataR, pNGradR, sizeR, nCompare, regWeight);

	free(pNGradR);
	free(pNGradL);

	return res;
}

#define MIN_SERIES_SIZE 20
RESULT_CODE myDynamicTimeWarpingPyramid(sPnt* pathDTW, int* sizePath, float* matDTW,
	int* pDataL, int sizeL,
	int* pDataR, int sizeR,
	int sizeSeriesWindow, int sizeMetricWindow, int* nScales)
{
	if (sizeL != sizeR) return ERROR_WRONG_INPUT;
	
	const int scaleStep = 2;
	int usedScales = 1, cSeriesSize = sizeL, scaleFactor = 1;

	for (; usedScales < *nScales; ++usedScales)
	{
		if (cSeriesSize < MIN_SERIES_SIZE)
			break;
		cSeriesSize /= 2; scaleFactor *= scaleStep;
	}
	std::cout << "LOG: Series size: " << cSeriesSize << ", scaleFactor: " << scaleFactor << ", usedScales: " << usedScales << std::endl;

	// float* matDTW = (float*)malloc(sizeL * sizeR * sizeof(float));
	int sizeP = sizeL + sizeR + 2;
	sPnt* path = (sPnt*)malloc(sizeP * sizeof(sPnt));
	int* pNDataR = (int*)malloc(sizeR * sizeof(int));
	int* pNGradR = (int*)malloc(sizeR * sizeof(int));
	int* pNDataL = (int*)malloc(sizeL * sizeof(int));
	int* pNGradL = (int*)malloc(sizeL * sizeof(int));
	
	for (int i = 0; i < usedScales; ++i)
	{
		std::cout << "LOG: Series size: " << cSeriesSize << ", scaleFactor: " << scaleFactor << std::endl;

		downsampleSeries(pNDataR, cSeriesSize, pDataR, sizeR, scaleFactor);
		sequenceGradient(pNGradR, pNDataR, cSeriesSize);
		
		downsampleSeries(pNDataL, cSeriesSize, pDataL, sizeL, scaleFactor);
		sequenceGradient(pNGradL, pNDataL, cSeriesSize);
		float* _matDTW = (float*)malloc(cSeriesSize * cSeriesSize * sizeof(float));

		int res;
		if (i == 0) // minimal level
		{
			std::cout << "LOG: Zero level: calculating full dtw" << std::endl;
			res = calcFullDynamicTimeWarping(path, &sizeP, _matDTW, pNDataL, pNGradL, cSeriesSize, pNDataR, pNGradR, cSeriesSize, -1, 0);
			for (int i = 0; i < sizeP; ++i) {
				std::cout << "(" << path[i].x << "," << path[i].y << ") ";
			}
			std::cout << std::endl;
			// break;
		}
		else // upper levels
		{
			std::cout << "LOG: level: " << i << " calculating local dtw" << std::endl;
			std::cout << "LOG: upsample res" << std::endl;
			upsamplePath(path, &sizeP, 2);
			res = calcFullDynamicTimeWarping(path, &sizeP, _matDTW, pNDataL, pNGradL, cSeriesSize, pNDataR, pNGradR, cSeriesSize, sizeSeriesWindow, 0);
		}
		if (i == usedScales - 1)
		{
			for (int k = 0; k < sizeR*sizeL; ++k)
				matDTW[k] = _matDTW[k];
		}

		free(_matDTW);

		if (ERROR_OK != res)
		{
			free(pNGradL);
			free(pNGradR);
			free(pNDataL);
			free(pNDataR);
			free(path);
			return res;
		}

		scaleFactor /= scaleStep;
		cSeriesSize *= scaleStep;
	}

	for (int i = 0; i < sizeP; ++i)
	{
		pathDTW[i].x = path[i].x; pathDTW[i].y = path[i].y;
	}
	free(pNGradL);
	free(pNGradR);
	free(pNDataL);
	free(pNDataR);
	free(path);
	return ERROR_OK;
}


RESULT_CODE myDynamicTimeWarping(sPnt* pathDTW, float* matDTW, int* sizePath, int* pDataL, int sizeL, int* pDataR, int sizeR, int sizeSeriesWindow, int sizeMetricWindow)
{
	int w_2 = sizeMetricWindow / 2; // metric window size
	short cw_2 = sizeSeriesWindow / 2;  // series comparison window size
	short sizePadL = sizeL + w_2, sizePadR = sizeR + w_2;
	
	// zero series padding
	int* pDataPadL = (int*)malloc(sizePadL * sizeof(int));
	int* pDataPadR = (int*)malloc(sizePadR * sizeof(int));
	
	memset(pDataPadL, 0, sizePadL * sizeof(int));
	memset(pDataPadR, 0, sizePadR * sizeof(int));
	
	memcpy(pDataPadL, pDataL, sizeL * sizeof(int));
	memcpy(pDataPadR, pDataR, sizeR * sizeof(int));

	// DTW matrix
	sPnt* matrixPrevState = (sPnt*)malloc(sizeL * sizeR * sizeof(sPnt));
	// infinity on edges
	for (int i = 0; i < sizeL; ++i)
	{
		for (int j = 0; j < sizeR; ++j)
		{
			matDTW[i * sizeR + j] = INFINITY;
			matrixPrevState[i * sizeR + j].x = (i < j) ? j - 1: j;
			matrixPrevState[i * sizeR + j].y = (i < j) ? i : i - 1;
		}
		if (i < w_2)
		{
			matDTW[i * sizeR + i] = MIA_abs(pDataR[i] - pDataPadL[i]);
			matrixPrevState[i * sizeR + i].x = i - 1;
			matrixPrevState[i * sizeR + i].y = i - 1;
		}
	}

	for (int i = w_2; i < sizeL; ++i)
	{
		
		int jBeg = MIA_max(w_2, i - cw_2), jEnd = MIA_min(sizeR-1, i + cw_2);
		float* ptrC = matDTW + (i)* sizeR + jBeg;
		float* ptrT = matDTW + (i - 1) * sizeR + jBeg;
		float* ptrL = matDTW + i * sizeR + jBeg - 1;
		float* ptrTL = matDTW + (i - 1) * sizeR + jBeg - 1;
		sPnt* ptrPrev = matrixPrevState + i * sizeR + jBeg;

		for (int j = jBeg; j <= jEnd; ++j)
		{
			*ptrC = metricWindowDTW(pDataPadL, i, pDataPadR, j, w_2);
			
			if (*ptrT <= MIA_min(*ptrL, *ptrTL))
			{
				*ptrC += *ptrT;
				ptrPrev->x = j;
				ptrPrev->y = i-1;
			}
			else if (*ptrL <= MIA_min(*ptrT, *ptrTL))
			{
				*ptrC += *ptrL;
				ptrPrev->x = j - 1;
				ptrPrev->y = i;
			}
			else
			{
				*ptrC += *ptrTL;
				ptrPrev->x = j-1;
				ptrPrev->y = i-1;
			}
			++ptrPrev; ++ptrC; ++ptrL; ++ptrT; ++ptrTL;
		}
	}
	free(pDataPadL);
	free(pDataPadR);

	// search for the minimal total cost
	sPnt optPathEnd = { -1, -1 };
	float optPathValue = INFINITY;
	float* ptrLastCol = matDTW + sizeR * sizeL - 1;
	float* ptrLastRow = ptrLastCol;
	// fprintf(stdout, "sizeL: %d, sizeR: %d\n", sizeL, sizeR);
	for (int i = 0; i < cw_2; ++i)
	{
		if (*ptrLastRow < optPathValue)
		{
			optPathValue = *ptrLastRow;
			optPathEnd.x = sizeR - i - 1;
			optPathEnd.y = sizeL - 1;
		}
		if (*ptrLastCol < optPathValue)
		{
			optPathValue = *ptrLastCol;
			optPathEnd.x = sizeR - 1;
			optPathEnd.y = sizeL - i - 1;
		}
		--ptrLastRow; ptrLastCol -= sizeR;
	}
	
	// obtain minimal total cost path
	for (int i = 0; i < *sizePath; ++i)
	{
		int xEnd = optPathEnd.x, yEnd = optPathEnd.y;
		// std::cout << "i: " << i << "xEnd: " << xEnd << " yEnd: " << yEnd << std::endl;
		if (xEnd <= 0 || yEnd <= 0)
		{
			*sizePath = i;
			break;
		}
		pathDTW[i].x = optPathEnd.x;
		pathDTW[i].y = optPathEnd.y;
		optPathEnd.x = matrixPrevState[yEnd * sizeR + xEnd].x;
		optPathEnd.y = matrixPrevState[yEnd * sizeR + xEnd].y;
	}
	
	for (int i = 0; i < *sizePath / 2; ++i)
	{
		sPnt t;
		t.x = pathDTW[i].x;
		t.y = pathDTW[i].y;
		pathDTW[i].x = pathDTW[*sizePath - i - 1].x;
		pathDTW[i].y = pathDTW[*sizePath - i - 1].y;
		pathDTW[*sizePath - i - 1].x = t.x;
		pathDTW[*sizePath - i - 1].y = t.y;
	}
	free(matrixPrevState);
	return ERROR_OK;
}

RESULT_CODE meanProjectionDynamicTimeWarping(int *pProjAvg, int* pProjL, int projSizeL, int* pProjR, int projSizeR)
{
	float* matDTW = (float*)malloc(projSizeL * projSizeR * sizeof(float));
	int sizePath = projSizeL + projSizeR + 2;
	sPnt *pathDTW = (sPnt*)malloc(sizePath * sizeof(sPnt));
	int sizeMetricWin = 3, sizeCompareWin = 10;
	int res = myDynamicTimeWarping(pathDTW, matDTW, &sizePath, pProjL, projSizeL, pProjR, projSizeR, sizeCompareWin, sizeMetricWin);
	if (ERROR_OK != res)
	{	
		free(matDTW);
		free(pathDTW);
		return res;
	}
	
	int projSize = 0;
	for (int i = 0; i < sizePath; ++i)
	{
		pProjAvg[pathDTW[i].y] = pProjL[pathDTW[i].y] + pProjR[pathDTW[i].x];
		if (++projSize == projSizeL)
			break;
	}

	free(matDTW);
	free(pathDTW);
	return 0;
}

RESULT_CODE myDetectBaseRadii(
    SPupilInfo* psPI,         // OUT: pupil data
    SIrisInfo* psII,          // OUT: iris data
    SCenterInfo* psCI,  // IN:  center data
    unsigned char* im,  // IN:  source image
    int W,                   // IN:  image size
    int H,                   //
    int* kernel3x3,           // edge detection kernel: kernel3x3[] = { -3, 0, 3, -10, 0, 10, -3, 0, 3 };
    int angMin,               // minimal projection estimation angle
    int angMax,               // maximal projection estimation angle
    int minR,                 // minimal radius
    int maxR,                 // maximal radius
    int thrHoriz,             // threshold for horizontal derivative
    float thrDot,             // threshold for circular shape: 0.8 - 1.0
    int* pProjR,              // circular projection - right side nums
    int* pProjL,               // circular projection - left side nums
	int blurHW
)
{
    int res = ERROR_OK;
    // int maxR = MIA_min(MIA_max(MIA_min(psCI->xc, W - psCI->xc), MIA_min(psCI->yc, H - psCI->yc)), H / 2);

    if ((res = FindHoughProjection(pProjR, pProjL, angMin, angMax, kernel3x3, thrHoriz,thrDot,minR,maxR, im, W, H, psCI->xc, psCI->yc, blurHW)))
    {
        fprintf(stderr, "[ERROR]: Projection feature generation failed.\n");
        return res;
    }
    memset(pProjL, 0, minR * sizeof(int));
    memset(pProjR, 0, minR * sizeof(int));

	for (int i = maxR; i < W; ++i)
	{
		pProjR[i] = 0; pProjL[i] = 0;
	}

    int mode = IVIR_ARRPOXIR_PRPRPR | IVIR_ARRPOXIR_USEPROJ | IVIR_ARRPOXIR_PRPRPR; //  IVIR_ARRPOXIR_LEFEYE | IVIR_ARRPOXIR_PRPRPR;
    if ((res = IVIR_PrPu(
        psPI, psII, psCI, im, W, H,
        mode, -1, -1, pProjR, pProjL, NULL)))
    {
        fprintf(stderr, "[ERROR]: base radii estimations failed.\n");
    }
    return res;
}


int DetectBaseRadii(
    SPupilInfo* psPI,         // OUT: pupil data
    SIrisInfo* psII,          // OUT: iris data
    const SCenterInfo* psCI,  // IN:  center data
    const unsigned char* im,  // IN:  source image
    int xs,                   // IN:  image size
    int ys,                   //
    int* pProjR,               // circular projection - right side nums
    int* pProjL               // circular projection - left side nums
                              // int mode // mask: 1- draw projs to buf, 2- use projs from buf
)
{
    RESULT_CODE res;
    int outRadius = xs / 2;
    int buflen = 0;

    if ((res = IPL_PROJ_FindHoughDonatorProjection7(
        pProjR, pProjL,
        im, xs, ys, psCI->xc, psCI->yc,
        MINIMAL_PUPIL_RADIUS, &outRadius, 0,// dense analysis
        NULL, &buflen)))
    {
        fprintf(stderr, "[ERROR]: buffer size calculation failed.\n");
        return res;
    }
    void* buf = malloc(buflen);
    if ((res = IPL_PROJ_FindHoughDonatorProjection7(
        pProjR, pProjL,
        im, xs, ys, psCI->xc, psCI->yc,
        MINIMAL_PUPIL_RADIUS, &outRadius, 0,// dense analysis
        buf, &buflen)))
    {
        free(buf);
        fprintf(stderr, "[ERROR]: projection computation fails.\n");
        return res;
    }

    int mode = IVIR_ARRPOXIR_PRPRPR; //  IVIR_ARRPOXIR_LEFEYE | IVIR_ARRPOXIR_PRPRPR;
    if ((res = IVIR_PrPu(
        psPI, psII, psCI, im, xs, ys,
        mode, -1, -1, pProjR, pProjL, NULL)))
    {
        free(buf);
        fprintf(stderr, "[ERROR]: base radii estimations failed.\n");
        return res;
    }
    free(buf);
    return ERROR_OK;
}
