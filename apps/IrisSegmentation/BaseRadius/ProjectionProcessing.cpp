/*----------------------------------------------------------------------------*/
/*                                                                            */
/*      Created:  5 May 2011                                                  */
/*      Revision: 1.0.00                                                      */
/*      Purpose:  bazrad detection by projections method                      */
/*      Authors:                                                              */
/*        Ivan Matveev   , Iurii Efimov										  */
/*                                                                            */
/*----------------------------------------------------------------------------*/
#define __FILENUM__ 198 // __FILENUM__TAG198

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <memory.h>
#include <iostream>

#include "GlobalParams.h"
#include "BaseRadii.h"
#include "ImageProcessing.h"
#include "BorderDetection.h"

typedef struct
{
	int idx;
	int pos;
	int val;
} SProjectionMaximum;

typedef struct
{
    int q;
    int xpl;
    int xpr;
    int xil;
    int xir;
    int vpl;
    int vpr;
    int vil;
    int vir;
    char casechar;
} SBazRadHyp;



RESULT_CODE FindHoughProjection(
    int* pProjR,
    int* pProjL,
    int angBeg,
    int angEnd,
    int* Kernel3x3,
    int thrHoriz,
    int thrDot,
    int minR,
    int maxR,
    unsigned char* im,
    int xs,
    int ys,
    int xc,
    int yc)
{
    if (angEnd > 90 || angEnd < -90 || angBeg < -90 || angBeg > 90 || angEnd <= angBeg)
    {
        fprintf(stderr, "ERROR: FindHoughProjection: wrong angle parameters.\n");
        return ERROR_WRONG_INPUT;
    }
    if (NULL == im || xs < 1 || ys < 1)
    {
        fprintf(stderr, "ERROR: FindHoughProjection: wrong src image parameters.\n");
        return ERROR_WRONG_INPUT;
    }
    if (NULL == Kernel3x3)
    {
        fprintf(stderr, "ERROR: FindHoughProjection: wrong edge kernel parameters.\n");
        return ERROR_WRONG_INPUT;
    }
    if (NULL == pProjL)
    {
        fprintf(stderr, "ERROR: FindHoughProjection: left projection buffer is NULL.\n");
        return ERROR_WRONG_INPUT;
    }
    if (NULL == pProjR)
    {
        fprintf(stderr, "ERROR: FindHoughProjection: right projection buffer is NULL.\n");
        return ERROR_WRONG_INPUT;
    }

    float kEnd = tanf((float)(angEnd) / 180 * PI), kBeg = tanf((float)(angBeg) / 180 * PI);
    int minY, maxY;
    // numerical limits
    if (kBeg <= 1e-4 && kBeg >= -1e-4)
        minY = yc;
    else
    {
        float xIntercept = -(float)yc / (kBeg + 1e-7) + xc;
        // printf("intercept: %f\n", xIntercept);
        minY = (int)(((xIntercept < (float)xs) ? .0 : -kBeg * (xIntercept - (float)xs)) + .5f);
    }

    if (kEnd <= 1e-4 && kEnd >= -1e-4)
        maxY = yc;
    else
    {
        float xIntercept = (float)(ys - yc) / (kEnd + 1e-7) + xc;
        // printf("intercept: %f\n", xIntercept);
        maxY = (int)(((xIntercept < (float)xs) ? (float)ys : (float)ys - kEnd * (xIntercept - (float)xs)) + .5f);
    }
    // sanity check
    minY = MIA_max(0, minY);
    maxY = MIA_min(maxY, ys);

    // num points in circle
    int* numPoR = (int*)malloc(xs * sizeof(int));
    int* numPoL = (int*)malloc(xs * sizeof(int));
    memset(numPoL, 0, xs * sizeof(int));
    memset(numPoR, 0, xs * sizeof(int));
    int* numPoR_b = (int*)malloc(xs * sizeof(int));
    int* numPoL_b = (int*)malloc(xs * sizeof(int));

    int* pnProjR = (int*)malloc(xs * sizeof(int));
    int* pnProjL = (int*)malloc(xs * sizeof(int));
    memset(pnProjR, 0, xs * sizeof(int));
    memset(pnProjL, 0, xs * sizeof(int));

    // right side
    for (int y = minY + 1; y < maxY - 1; ++y)
    {
        uint8 *pImRow = im + y * xs + xc - 1;
        uint8 *pImRowPrev = im + (y - 1) * xs + xc - 1;
        uint8 *pImRowNext = im + (y + 1) * xs + xc - 1;

        for (int x = 0; x < xs - xc - 1; ++x)
        {
            ++pImRow;
            ++pImRowPrev;
            ++pImRowNext;

            int r = INTSQRT(x*x + (y - yc)*(y - yc));

            int gx =
                Kernel3x3[0] * pImRowPrev[-1] + Kernel3x3[1] * pImRowPrev[0] + Kernel3x3[2] * pImRowPrev[1] +
                Kernel3x3[3] * pImRow[-1] + Kernel3x3[4] * pImRow[0] + Kernel3x3[5] * pImRow[1] +
                Kernel3x3[6] * pImRowNext[-1] + Kernel3x3[7] * pImRowNext[0] + Kernel3x3[8] * pImRowNext[1];

            // std ::cout << gx << " ";
            if (gx < thrHoriz)
                continue;
            ++numPoR[r];
            if ((float)x * kBeg <= y - yc && (float)x * kEnd >= y - yc)
            {
                // *pImRow = 255;
                int gy =
                    Kernel3x3[0] * pImRowPrev[-1] + Kernel3x3[3] * pImRowPrev[0] + Kernel3x3[6] * pImRowPrev[1] +
                    Kernel3x3[1] * pImRow[-1] + Kernel3x3[4] * pImRow[0] + Kernel3x3[7] * pImRow[1] +
                    Kernel3x3[2] * pImRowNext[-1] + Kernel3x3[5] * pImRowNext[0] + Kernel3x3[8] * pImRowNext[1];
                // std::cout << gy << " ";
                int G = INTSQRT(gx * gx + gy * gy);
                // std::cout << G << " " << r;
                if ((float)(gx * x + gy * (y - yc)) >= 2 * thrDot * G * r)
                {
                    // std::cout << gx << " " << gy << " " << G  << " " << r << std:: endl;
                    ++pnProjR[r];
                }
            }
        }
        // std::cout << std::endl;
    }
    // left side
    kBeg = -kBeg;
    kEnd = -kEnd;
    // std::cout << "Left side!" << kBeg << " " << kEnd << std::endl;
    if (kBeg <= 1e-4 && kBeg >= -1e-4)
        minY = yc - 1;
    else
    {
        float xIntercept = -(float)yc / (kBeg + 1e-7) + xc;
        // printf("intercept: %f\n", xIntercept);
        minY = (int)(((xIntercept >= 0) ? .0 : -kBeg * xIntercept + .5f));
    }

    if (kEnd <= 1e-4 && kEnd >= -1e-4)
        maxY = yc + 1;
    else
    {
        float xIntercept = (float)(ys - yc) / (kEnd + 1e-7) + xc;
        // printf("intercept: %f\n", xIntercept);
        maxY = (int)(((xIntercept >= 0) ? (float)ys : (float)ys - kEnd * xIntercept + .5f));
    }
    // sanity check
    minY = MIA_max(0, minY);
    maxY = MIA_min(maxY, ys);
    // std::cout << "Min Y " << minY << "Max Y " << maxY << std::endl;
    for (int y = minY + 1; y < maxY - 1; ++y)
    {
        uint8 *pImRow = im + y * xs + xc + 1;
        uint8 *pImRowPrev = im + (y - 1) * xs + xc + 1;
        uint8 *pImRowNext = im + (y + 1) * xs + xc + 1;

        int xBeg = (y < yc) ? (int)((float)(y - yc) / kBeg + 1) : (int)((float)(y - yc) / kEnd + 1);
        for (int x = 0; x > -xc; --x)
        {
            --pImRow;
            --pImRowPrev;
            --pImRowNext;

            int r = INTSQRT(x*x + (y - yc)*(y - yc));

            int gx =
                Kernel3x3[0] * pImRowPrev[-1] + Kernel3x3[1] * pImRowPrev[0] + Kernel3x3[2] * pImRowPrev[1] +
                Kernel3x3[3] * pImRow[-1] + Kernel3x3[4] * pImRow[0] + Kernel3x3[5] * pImRow[1] +
                Kernel3x3[6] * pImRowNext[-1] + Kernel3x3[7] * pImRowNext[0] + Kernel3x3[8] * pImRowNext[1];
            if (gx > -thrHoriz)
                continue;

            ++numPoL[r];

            if ((float)x * kBeg <= (float)(y - yc) && (float)x * kEnd >= (float)(y - yc))
            {
                // *pImRow = 254;
                int gy =
                    Kernel3x3[0] * pImRowPrev[-1] + Kernel3x3[3] * pImRowPrev[0] + Kernel3x3[6] * pImRowPrev[1] +
                    Kernel3x3[1] * pImRow[-1] + Kernel3x3[4] * pImRow[0] + Kernel3x3[7] * pImRow[1] +
                    Kernel3x3[2] * pImRowNext[-1] + Kernel3x3[5] * pImRowNext[0] + Kernel3x3[8] * pImRowNext[1];

                int G = INTSQRT(gx * gx + gy * gy);
                if (gx * x + gy * (y - yc) >= 2 * thrDot * G * r)
                {
                    // std::cout << r << " " << pnProjL[r] << std::endl;
                    ++pnProjL[r];
                }
            }
        }
    }

    // normalize
    for (int i = minR; i <= maxR; i++)
    {
        // std::cout << numPoR[i] << " ";
        // pnProjR[i] = (numPoR[i]) ? (int)((1024 * pnProjR[i]) / numPoR[i]) : 0;
        pnProjR[i] = (numPoR[i]) ? (int)((1024 * pnProjR[i]) / (i)) : 0;
        // std::cout << (float)(pnProjR[i]) / (numPoR[i] + 1e-5) << std::endl;
        pnProjL[i] = (numPoL[i]) ? (int)((1024 * pnProjL[i]) / (i)) : 0;
    }
    for (int i = 0; i < minR; i++)
    {
        numPoR[i] = numPoR[minR];
        numPoL[i] = numPoL[minR];
        pnProjR[i] = pnProjR[minR];
        pnProjL[i] = pnProjL[minR];
    }
    // std::cout << std::endl;
    // process histogram: blur
    IPL_HIST_Blur(numPoR_b, numPoR, maxR + 1, BLUR_HIST_WND);
    IPL_HIST_Blur(numPoL_b, numPoL, maxR + 1, BLUR_HIST_WND);
    IPL_HIST_Blur(pProjR, pnProjR, maxR + 1, BLUR_HIST_WND);
    IPL_HIST_Blur(pProjL, pnProjL, maxR + 1, BLUR_HIST_WND);

    // correct for small presence
    for (int i = 1; i <= maxR; i++)
    {
        int L = 10000 * numPoR_b[i] / (157 * (2 * BLUR_HIST_WND + 1)*i);
        if (L < 20)
            pProjR[i] = 0;
        L = 10000 * numPoL_b[i] / (157 * (2 * BLUR_HIST_WND + 1)*i);
        if (L < 20)
            pProjL[i] = 0;
    }

    free(pnProjL);
    free(pnProjR);
    free(numPoR);
    free(numPoL);
    free(numPoR_b);
    free(numPoL_b);
    return ERROR_OK;
}



static int ownIVIR_SortMaxima_by_ValPos(
	const void* e1,
	const void* e2)
{
	if (((SProjectionMaximum*)e2)->val - ((SProjectionMaximum*)e1)->val)
		return ((SProjectionMaximum*)e2)->val - ((SProjectionMaximum*)e1)->val;
	return ((SProjectionMaximum*)e2)->pos - ((SProjectionMaximum*)e1)->pos;
}

static int ownIVIR_SortMaxima_by_Pos(
	const void* e1,
	const void* e2)
{
	return ((SProjectionMaximum*)e1)->pos - ((SProjectionMaximum*)e2)->pos;
}

void IVIR_ProjPreProcess(
	int* hist,
	int len)
{
	int *pro2 = NULL, i, j, Sum, k, Num;

	pro2 = (int*)malloc(len * sizeof(pro2[0]));
	for (i = 0; i<len; i++)
	{
		Sum = Num = 0;
		for (j = -len / 8; j <= len / 8; j++)
		{
			k = i + j;
			if (k<0)
				continue;
			if (k >= len)
				continue;
			Num++;
			Sum += hist[k];
		}
		pro2[i] = (Num) ? (Sum / (Num * 2)) : 0;
	}
	for (i = 0; i<len; i++)
		hist[i] = (hist[i]>pro2[i]) ? (hist[i] - pro2[i]) : 0;
	free(pro2);
}

static void IVIR_ProjPostProcess2(
	int* hist,
	int len,
	int norm)
{
	int i;

	if (norm <= 0)
		return;
	for (i = 0; i<len; i++)
		if (hist[i] <= norm)
			hist[i] = (hist[i] * 1000) / norm;
		else
			hist[i] = 999;
}

// based on pupil info, find projection local maxima positions suitable for iris
/*RESULT_CODE IVIR_FindMaximaPositions(
	int* pMaxes,        // OUT: positions of local maxima (sorted by value)
	int* vMaxes,        // OUT: values of local maxima
	int** pProjL,       // OUT: pointer to left projection buffer (out of pvbuf)
	int** pProjR,       // OUT: pointer to right projection buffer (out of pvbuf)
	int* pnProjLNorm,   // OUT: left projection norm
	int* pnProjRNorm,   // OUT: right projection norm
	int r,              // IN:  inner radius
	int* pR,            // IN/OUT: outer radius
	const unsigned char* im,  // IN:  source image
	int xs,                   // IN:  image size
	int ys,                   //
	const SPupilInfo* pPupil, // IN:  pupil data
	int nMaxes,               // IN:  number of required maxima
	void* pvBuf,              // IN:  temporary buffer
	int* pnLen,               // IN/OUT: number of bytes allocated/used
	const char* nambeg)       // IN:  for debugging purposes, set to NULL
{
	int nLoc, sz, i;
	RESULT_CODE res;
	int *projL, *projR, nMaxL, nMaxR, *anMaxsL, *anMaxsR, dummy_R;
	void* projbuf;
	SProjectionMaximum *LefMax, *RigMax;

	// check arguments
	if ((xs <= 0) || (ys <= 0) || (nMaxes <= 0))
		return ERROR_WRONG_INPUT;
	if (pPupil == NULL)
		return ERROR_NULL_POINTER;
	if ((pPupil->q1 <= 0) || (pPupil->r <= 0))
		return ERROR_WRONG_INPUT;
	// get maximum radii of ring
	if (r <= 0)
		r = pPupil->r * 4 / 3;
	if (pR == NULL)
	{
		pR = &dummy_R;
		*pR = pPupil->r * 6;
	}
	res = IPL_PROJ_FindHoughDonatorProjection6(
		NULL, NULL, NULL, xs, ys, pPupil->xc, pPupil->yc,
		r, pR, 0, NULL, &nLoc);
	if (res != ERROR_OK)
		return res;
	// add this function's own requirements
	sz = nLoc +                // donator temp mem
		2 * sizeof(int)*xs +    // projections
		2 * sizeof(int)*xs +      // maxima values
		sizeof(SProjectionMaximum)*xs; // maxima structures
									   // check size
	if (pvBuf == NULL)
	{
		*pnLen = sz;
		return ERROR_OK;
	}
	// check more args
	if ((pMaxes == NULL) || (vMaxes == NULL) || (im == NULL))
		return ERROR_NULL_POINTER;
	//    for (i=0;i<nMaxes*2;i++)
	//    pMaxes[i] = vMaxes[i] = 0;
	// allocate
	projR = (int*)pvBuf;    // projections
	projL = projR + xs;
	anMaxsR = projL + xs;    // maxima values
	anMaxsL = anMaxsR + xs;
	LefMax = (SProjectionMaximum*)(anMaxsL + xs); // maxima structures
	RigMax = LefMax + xs / 2;
	projbuf = (void*)(RigMax + xs / 2); // donator temp mem
	if (pProjL)
		*pProjL = projL;
	if (pProjR)
		*pProjR = projR;
	// calculate left side, right side and total circular projection in ring
	res = IPL_PROJ_FindHoughDonatorProjection6(
		projR,          // OUT: circular projection - right side nums
		projL,          // OUT: circular projection - left side nums
		im,  // IN:  image
		xs,                   // IN:  image size
		ys,                   // 
		pPupil->xc,                   // IN:  ring center
		pPupil->yc,                   //
		r,                    // IN:  inner radius
		pR,                  // IN/OUT: outer radius proposed/actually used
		0,       // IN:  maximum radius for dense processing. 0 - dense always
		projbuf,                // IN:  external buffer
		&nLoc);              // IN/OUT: allocated/used bytes
	if (res != ERROR_OK)
		return res;
	IVIR_ProjPreProcess(projR, *pR);
	IVIR_ProjPreProcess(projL, *pR);
	// find max values
	nMaxL = ownIVIR_ListMaxs(anMaxsL, projL, r, *pR - 1);
	nMaxR = ownIVIR_ListMaxs(anMaxsR, projR, r, *pR - 1);
	// fill in maxima structures
	for (i = 0; i<nMaxL; i++)
	{
		LefMax[i].idx = i;
		LefMax[i].pos = anMaxsL[2 * i];
		LefMax[i].val = projL[anMaxsL[2 * i]];
	}
	for (i = 0; i<nMaxR; i++)
	{
		RigMax[i].idx = i;
		RigMax[i].pos = anMaxsR[2 * i];
		RigMax[i].val = projR[anMaxsR[2 * i]];
	}
	// sort maxima in value_descending/coord_descending order
	qsort(LefMax, nMaxL, sizeof(LefMax[0]), ownIVIR_SortMaxima_by_ValPos);
	qsort(RigMax, nMaxR, sizeof(RigMax[0]), ownIVIR_SortMaxima_by_ValPos);
	// cut number of maxima
	if (nMaxL>nMaxes)
		nMaxL = nMaxes;
	if (nMaxR>nMaxes)
		nMaxR = nMaxes;
	// sort maxima in coord_ascending order
	//    qsort(LefMax,nMaxL,sizeof(LefMax[0]),ownIVIR_SortMaxima_by_Pos);
	//  qsort(RigMax,nMaxR,sizeof(RigMax[0]),ownIVIR_SortMaxima_by_Pos);
	// if number of located maxes is less than required, fill with negs
	for (i = nMaxL; i<nMaxes; i++)
		LefMax[i].idx = LefMax[i].pos = LefMax[i].val = -1;
	for (i = nMaxR; i<nMaxes; i++)
		RigMax[i].idx = RigMax[i].pos = RigMax[i].val = -1;
	if (pnProjLNorm)
		*pnProjLNorm = projL[LefMax[0].pos];
	if (pnProjRNorm)
		*pnProjRNorm = projR[RigMax[0].pos];
	IVIR_ProjPostProcess2(projL, *pR, projL[LefMax[0].pos]);
	IVIR_ProjPostProcess2(projR, *pR, projR[RigMax[0].pos]);
	// output
	for (i = 0; i<nMaxes; i++)
	{
		pMaxes[i] = LefMax[i].pos;
		vMaxes[i] = (LefMax[i].pos>0) ? (projL[LefMax[i].pos]) : 0;
		pMaxes[i + nMaxes] = RigMax[i].pos;
		vMaxes[i + nMaxes] = (RigMax[i].pos>0) ? (projR[RigMax[i].pos]) : 0;
	}
	if (nambeg)
	{
		unsigned char* drawim;
		char nama[FILENAME_MAX];
		drawim = (unsigned char*)malloc(xs*ys * 3);
		DBGL_PUMP_Make3Bpp(drawim, xs * 3, im, xs, xs, ys);
		// center, search bound, true pupil, true iris
		DBGL_DRAW_CircleInRGB(drawim, xs * 3, xs, ys, pPupil->xc, pPupil->yc, r, 0x00ff00);
		DBGL_DRAW_CircleInRGB(drawim, xs * 3, xs, ys, pPupil->xc, pPupil->yc, *pR, 0x00ff00);
		DBGL_DRAW_CircleInRGB(drawim, xs * 3, xs, ys, pPupil->xc, pPupil->yc, pPupil->r, 0x003f00);
		// histograms with maxima and principal maximum
		DBGL_DRAW_HistogramInRGBext(drawim, *pR, 100, xs * 3, projL, *pR, 1);
		DBGL_DRAW_HistogramInRGBext(drawim + (xs - *pR) * 3, *pR, 100, xs * 3, projR, *pR, 1);
		DBGL_DRAW_StrobesInRGB(drawim, *pR, 100, xs * 3, projL, *pR, anMaxsL, nMaxL * 2, 0xff0000);
		DBGL_DRAW_StrobesInRGB(drawim + (xs - *pR) * 3, *pR, 100, xs * 3, projR, *pR, anMaxsR, nMaxR * 2, 0xff0000);
		for (i = 0; i<nMaxL; i++)
			anMaxsL[i] = LefMax[i].pos;
		for (i = 0; i<nMaxR; i++)
			anMaxsR[i] = RigMax[i].pos;
		DBGL_DRAW_StrobesInRGB(drawim, *pR, 100, xs * 3, projL, *pR, anMaxsL, nMaxL, 0x0000ff);
		DBGL_DRAW_StrobesInRGB(drawim + (xs - *pR) * 3, *pR, 100, xs * 3, projR, *pR, anMaxsR, nMaxR, 0x0000ff);
		// text report
		strcpy(nama, "Lef:");
		for (i = 0; i<nMaxL; i++)
			sprintf(nama + strlen(nama), "%d ", LefMax[i].pos);
		DBGL_TYPE_TextInRGB(drawim, xs, ys, nama, 0, ys / 2, 0xff00ff, -1);
		strcpy(nama, "Rig:");
		for (i = 0; i<nMaxL; i++)
			sprintf(nama + strlen(nama), "%d ", RigMax[i].pos);
		DBGL_TYPE_TextInRGB(drawim, xs, ys, nama, 0, ys / 2 + 10, 0xff00ff, -1);
		// file
		sprintf(nama, "%s_FMP.bmp", nambeg);
		DBGL_FILE_SaveRGB8Image(drawim, xs, xs, ys, nama);
		free(drawim);
	}
    
	return res;
}*/


static int ownIVIR_CorrectBazradHyp(
    int isLeft,
    const SBazRadHyp* psbrh,
    int nHyp,
    int Xpl,
    int Xpr,
    int Xil,
    int Xir,
    const int* projR,
    const int* anMaxsR,
    int nMaxR,
    const int* projL,
    const int* anMaxsL,
    int nMaxL,
    int ProjLen)
{
    int hypidx, maxidxNew, IrNew, maxidxOld;//VolOld,VolNew,

    // find hyp
    for (hypidx = 0; hypidx<nHyp; hypidx++)
        if ((psbrh[hypidx].xpl == Xpl) && (psbrh[hypidx].xpr == Xpr) &&
            (psbrh[hypidx].xil == Xil) && (psbrh[hypidx].xir == Xir))
            break;
    if (hypidx == nHyp)
        return -1;  // no hyp found
    if (!isLeft)
    { // correction for right
      // locate max on the right of hyp
        for (maxidxOld = 0; (maxidxOld<nMaxR) && (anMaxsR[2 * maxidxOld]<Xir); maxidxOld++);
        if (maxidxOld == nMaxR)
            return -2;  // cannot find index for old maximum
        maxidxNew = maxidxOld + 1;
        if (maxidxNew == nMaxR)
            return -3;  // no max outside of iris found
        IrNew = anMaxsR[2 * maxidxNew];
        //      if (2*(IrNew-Xir)>Xir-Xpr)
        //      return 0;   // too big correction required
        //      if (projR[IrNew]>2*projR[anMaxsR[2*maxidxNew-1]])
        //      return 0;   // gap between maxima is too deep
        //      if (projR[IrNew]*10>projR[Xir]*9)
        if (projR[IrNew]>projR[Xir])
            return IrNew;
        /*      if (projR[IrNew]*5>projR[Xir]*4)
        { // enough big projection value - decide
        VolOld = ownIVIR_GetPeakVol(projR,anMaxsR,nMaxR,maxidxOld,ProjLen);
        VolNew = ownIVIR_GetPeakVol(projR,anMaxsR,nMaxR,maxidxNew,ProjLen);
        return (VolOld>=VolNew)?0:IrNew;
        }*/
        return 0;
    }
    else
    { // correction for right
      // locate max on the right of hyp
        for (maxidxOld = 0; (maxidxOld<nMaxL) && (anMaxsL[2 * maxidxOld]<Xil); maxidxOld++);
        if (maxidxOld == nMaxL)
            return -2;  // cannot find index for old maximum
        maxidxNew = maxidxOld + 1;
        if (maxidxNew == nMaxL)
            return -3;  // no max outside of iris found
        IrNew = anMaxsL[2 * maxidxNew];
        //      if (2*(IrNew-Xil)>Xil-Xpl)
        //      return 0;   // too big correction required
        //      if (projL[IrNew]>2*projL[anMaxsL[2*maxidxNew-1]])
        //      return 0;   // gap between maxima is too deep
        //      if (projL[IrNew]*10>projL[Xil]*9)
        if (projL[IrNew]>projL[Xil])
            return IrNew;
        /*      if (projL[IrNew]*5>projL[Xil]*4)
        { // enough big projection value - decide
        VolOld = ownIVIR_GetPeakVol(projL,anMaxsL,nMaxL,maxidxOld,ProjLen);
        VolNew = ownIVIR_GetPeakVol(projL,anMaxsL,nMaxL,maxidxNew,ProjLen);
        return (VolOld>=VolNew)?0:IrNew;
        }*/
        return 0;
    }
}



int IVIR_ProjPostProcess(  // returns unnormed histogram maximum
	int* hist,
	int len)
{
	int *pro2 = NULL, i, j, Sum, k, Num;

	pro2 = (int*)malloc(len * sizeof(pro2[0]));
	for (i = 0; i<len; i++)
	{
		Sum = Num = 0;
		for (j = -len / 8; j <= len / 8; j++)
		{
			k = i + j;
			if (k<0)
				continue;
			if (k >= len)
				continue;
			Num++;
			Sum += hist[k];
		}
		pro2[i] = (Num) ? (Sum / (Num * 2)) : 0;
	}
	for (i = 0; i<len; i++)
		hist[i] = (hist[i]>pro2[i]) ? (hist[i] - pro2[i]) : 0;
	free(pro2);
	for (i = j = 0; i<len; i++)
		if (j<hist[i])
			j = hist[i];
	if (j <= 0)
		return j;
	for (i = 0; i<len; i++)
		hist[i] = (hist[i] * 1000) / j;
	return j;
}

int IVIR_P2I_GeomQuality(
	int PupilR,
	int posL,
	int posR)
{
	int nDecentration, QualDec, nPupirrat, Qualrat, nShift, QualShift, Qual;

	// check geometrical validity
	// decentration
	nDecentration = 100 * (posR - posL) / (posR + posL);
	if (nDecentration<0)
		nDecentration = -nDecentration;
	if (nDecentration>15)
		return 0; // invalid
				  /*if (nDecentration<=5)
				  QualDec = 100;
				  else
				  if (nDecentration<=10)
				  QualDec = 75+(nDecentration-10)*(100-75)/(5-10);
				  else
				  QualDec = 25+(nDecentration-15)*(75-25)/(10-15);*/
	if (nDecentration <= 8)
		QualDec = 100;
	else
		if (nDecentration <= 12)
			QualDec = 80 + (nDecentration - 12)*(100 - 80) / (8 - 12);
		else
			QualDec = 30 + (nDecentration - 15)*(80 - 30) / (12 - 15);
	// shift (another bounding for decentration)
	if (posL>posR)
		nShift = 100 * (posR - PupilR) / (posL - PupilR);
	else
		nShift = 100 * (posL - PupilR) / (posR - PupilR);
	if (nShift<60)
		return 0; // invalid
	if (nShift<75)
		QualShift = (nShift - 50) * 4;
	else
		QualShift = 100;
	// pupil-iris ratio
	nPupirrat = 200 * PupilR / (posR + posL);
	if ((nPupirrat >= 75) || (nPupirrat <= 13))
		return 0;
	if (nPupirrat<20)
		Qualrat = 0 + (nPupirrat - 12)*(100 - 0) / (20 - 12);
	else
		if (nPupirrat<60)
			Qualrat = 100;
		else
			if (nPupirrat<65)
				Qualrat = 100 + (nPupirrat - 60)*(75 - 100) / (65 - 60);
			else
				if (nPupirrat<70)
					Qualrat = 75 + (nPupirrat - 65)*(37 - 75) / (70 - 65);
				else
					Qualrat = 37 + (nPupirrat - 70)*(0 - 37) / (75 - 70);
	// common quality
	Qual = QualDec;
	if (Qual>Qualrat)
		Qual = Qualrat;
	if (Qual>QualShift)
		Qual = QualShift;
	return Qual;
}

/*typedef struct
{
    int idx;    // [IURII] originally missed
    int pos;    // position
    int val;    // value
    int sig;    // sign of gradient direction
    int wid;    // width of a peak at half-height
} SProjectionMaximum;*/


#define BLUR_HIST_WND 4

static int ownIVIR_difil(
    int a,
    int b)
{
    a = 1000 * (a - b) / (a + b);
    return (a >= 0) ? a : (-a);
}

static int ownIVIR_CutMax(
    int* nPos,
    const int* hist,
    int MaxCnt,
    int cutidx)
{
    int k;

    cutidx = (hist[nPos[cutidx * 2 - 1]]>hist[nPos[cutidx * 2 + 1]]) ? (cutidx * 2 - 1) : (cutidx * 2);
    // cut it away
    for (k = cutidx; k <= MaxCnt * 2 - 1; k++)
        nPos[k] = nPos[k + 2];
    return MaxCnt - 1;
}

static int ownIVIR_ListMaxs(
    int* nPos,
    const int* hist1,
    int minrad,
    int len)
{
    int MaxCnt, i, PossibleMaxBeg, PossibleMinBeg, Maxval;

    // === form sequence of max-min-max-...-max
    MaxCnt = 0;
    PossibleMaxBeg = PossibleMinBeg = 0;
    for (i = 1; i<len; i++)
    {
        if (hist1[i]>hist1[i - 1])
        { // growing
          // may be a new maximum at the point
            PossibleMaxBeg = i;
            // check if was minimum at previous point
            if ((PossibleMinBeg >= 0) && MaxCnt)
                nPos[MaxCnt * 2 - 1] = (PossibleMinBeg + i - 1) / 2;
            PossibleMinBeg = -1;
        }
        else
            if (hist1[i]<hist1[i - 1])
            { // falling
              // may be a new mimimum at the point
                PossibleMinBeg = i;
                // check if was maximum at previous point
                if (PossibleMaxBeg >= 0)
                {
                    nPos[MaxCnt * 2] = (PossibleMaxBeg + i - 1) / 2;
                    MaxCnt++;
                }
                PossibleMaxBeg = -1;
            }
    }
    if (PossibleMaxBeg >= 0)
    {
        nPos[MaxCnt * 2] = (PossibleMaxBeg + i - 1) / 2;
        MaxCnt++;
    }
    // minimum from the left
    nPos[-1] = 0;
    for (i = 1; i<nPos[0]; i++)
        if (hist1[nPos[-1]]>hist1[i])
            nPos[-1] = i;
    // minimum from the right
    nPos[2 * MaxCnt - 1] = len - 1;
    for (i = nPos[2 * (MaxCnt - 1)] + 1; i<len - 2; i++)
        if (hist1[nPos[2 * MaxCnt - 1]]>hist1[i])
            nPos[2 * MaxCnt - 1] = i;
    //== find maxval
    Maxval = hist1[nPos[0]];
    for (i = 1; i<MaxCnt; i++)
        if (Maxval<hist1[nPos[2 * i]])
            Maxval = hist1[nPos[2 * i]];
    // === sift maxima
    // single-max operations
    for (i = 0; i<MaxCnt; i++)
        if ((nPos[i * 2]<minrad) || (nPos[i * 2]>len - BLUR_HIST_WND))
        { // maximum position out of allowed range
            MaxCnt = ownIVIR_CutMax(nPos, hist1, MaxCnt, i);
            i--;
            continue;
        }
    // double-max operations
    for (i = 0; i<MaxCnt - 1; i++)
        if (nPos[(i + 1) * 2] - nPos[i * 2]<2 * BLUR_HIST_WND + 1)
        { // too close maximum
            MaxCnt = ownIVIR_CutMax(nPos, hist1, MaxCnt,
                (hist1[nPos[(i + 1) * 2]] <= hist1[nPos[i * 2]]) ? (i + 1) : i);
            i--;
            continue;
        }
    return MaxCnt;
}

void ownIVIR_GetMaxSubMax(
    int *pMax1,         // OUT: index of first max_submax
    int *pMax2,         // OUT: index of second max_submax
    const int* projL,
    const int* anMaxsL,
    int nMaxL,
    int *pMax3)
{
    int _max, _submax, i, _sbbmax;

    _max = 0;
    for (i = 1; i<nMaxL; i++)
        if (projL[anMaxsL[2 * _max]]<projL[anMaxsL[2 * i]])
            _max = i;
    _submax = -1;
    for (i = 0; i<nMaxL; i++)
        if ((anMaxsL[2 * i] * 4<anMaxsL[2 * _max] * 3) ||
            (anMaxsL[2 * i] * 3>anMaxsL[2 * _max] * 4))
            //      if ( ((anMaxsL[2*i]*4<anMaxsL[2*_max]*3)||
            //          (anMaxsL[2*i]*3>anMaxsL[2*_max]*4)  ) &&
            //       (projL[anMaxsL[2*_max]]<2*projL[anMaxsL[2*i]]) )
        {
            if (_submax == -1)
                _submax = i;
            else
                if (projL[anMaxsL[2 * _submax]]<projL[anMaxsL[2 * i]])
                    _submax = i;
        }
    if (_max<_submax)
    {
        *pMax1 = _max;
        *pMax2 = _submax;
    }
    else
    {
        *pMax1 = _submax;
        *pMax2 = _max;
    }
    *pMax3 = _sbbmax = -1;
    if (_submax == -1)
        return;
    for (i = 0; i<nMaxL; i++)
        if (((anMaxsL[2 * i] * 4<anMaxsL[2 * _max] * 3) || (anMaxsL[2 * i] * 3>anMaxsL[2 * _max] * 4)) &&
            ((anMaxsL[2 * i] * 4<anMaxsL[2 * _submax] * 3) || (anMaxsL[2 * i] * 3>anMaxsL[2 * _submax] * 4)))
        {
            if (_sbbmax == -1)
                _sbbmax = i;
            else
                if (projL[anMaxsL[2 * _sbbmax]]<projL[anMaxsL[2 * i]])
                    _sbbmax = i;
        }
    *pMax3 = _sbbmax;
}

int ownIVIR_sortSBRH(
    const void* e1,
    const void* e2)
{
    SBazRadHyp* h1 = (SBazRadHyp*)e1;
    SBazRadHyp* h2 = (SBazRadHyp*)e2;

    if (h2->q - h1->q)
        return h2->q - h1->q;
    if (h1->xpl - h2->xpl)
        return h1->xpl - h2->xpl;
    if (h1->xpr - h2->xpr)
        return h1->xpr - h2->xpr;
    if (h1->xil - h2->xil)
        return h1->xil - h2->xil;
    return h1->xir - h2->xir;
}

// detect approximate pupil and iris by projection method
#define __FILENUM__ 9 // __FILENUM__TAG9
#define BLUR_HIST_WND 4


// with blurring
RESULT_CODE IPL_PROJ_FindHoughDonatorProjection7(
  int* pnProjRN_b,          // OUT: circular projection - right side nums
  int* pnProjLN_b,          // OUT: circular projection - left side nums
  const unsigned char* im,  // IN:  image
  int xs,                   // IN:  image size
  int ys,                   // 
  int xc,                   // IN:  ring center
  int yc,                   //
  int r,                    // IN:  inner radius
  int* pR,                  // IN/OUT: outer radius proposed/actually used
  int nMaxRadOfDence,       // IN:  maximum radius for dense processing. 0 - dense always
  void* buf,                // IN:  external buffer
  int* buflen)              // IN/OUT: allocated/used bytes
{
  int i,j,gx,gy,y,G,L,S,*numpoR,*numpoL,currad,R;
  int j_inn,j_out,Y,D,RR_ii,rr_ii;
  int *pnProjRN,*pnProjLN,*numpoR_b,*numpoL_b;
  short *smoothedim,*psim;
  double* tempim;

    // check arguments
    if ((buflen==NULL)||(pR==NULL))
      return ERROR_NULL_POINTER;
    // verify radii and calculate max outer radius
    R = *pR;
    if (R<10)
      return ERROR_WRONG_INPUT;

    if ((r+5>R)||(r<5))
      return ERROR_WRONG_INPUT;

    if ((xs<2*r)||(ys<2*r))
      return ERROR_WRONG_INPUT;

    if ((xc-r/2<0)||(xc+r/2>=xs)||(yc-r/2<0)||(yc+r/2>=ys))
      return ERROR_WRONG_INPUT;
    // verify outer radius - bound by maximum span from center to border
    i = xs-1-xc;
    if (i<xc)
      i = xc;
    if (R>i)
      R = i;
    // verify outer radius - bound by image size
    if (2*R>xs)
      R = xs/2;
    // return corrected outer radius
    *pR = R;
    // calc size
    i = (2+4+2)*(R+1)*sizeof(numpoR[0])+
        xs*ys*sizeof(smoothedim[0])+xs*ys*sizeof(tempim[0]);
    if (buf==NULL)
    {
      *buflen = i;
      return ERROR_OK;
    }
    if (*buflen<i)
    {
      *buflen = i;
      return ERROR_MEMORY;// ERR_GEN_NOMEMORY;
    }
    // more check arguments
    if ((pnProjRN_b==NULL)||(pnProjLN_b==NULL)||(im==NULL))
      return ERROR_NULL_POINTER;
    /*if ((((ptr_t)im)&3)||(((ptr_t)buf)&3)||(xs&3)||
        (((ptr_t)pnProjRN_b)&3)||(((ptr_t)pnProjLN_b)&3))
      return ERROR_NULL_POINTER;
    */
    // allocate
    tempim = (double*)buf;
    numpoR   = (int*)(tempim+xs*ys);
    numpoL   = numpoR+(R+1);
    pnProjRN = numpoL+(R+1);
    pnProjLN = pnProjRN+(R+1);
    numpoR_b = pnProjLN+(R+1);
    numpoL_b = numpoR_b+(R+1);
    smoothedim = (short*)(numpoL_b+(R+1));
    // clear junk
    memset(pnProjRN,0,sizeof(pnProjRN[0])*(1+R));
    memset(pnProjLN,0,sizeof(pnProjLN[0])*(1+R));
    memset(numpoR,0,sizeof(numpoR[0])*(1+R));
    memset(numpoL,0,sizeof(numpoL[0])*(1+R));
    IPL_FILT_GaussianSmooth_uint8(
      smoothedim,
      tempim,
      im,
      xs,
      ys,
      2.);

    // along y
    if (!nMaxRadOfDence)
      D = 0x10000;
    else
    {
      D = (r*0x10000)/nMaxRadOfDence;
      if (D<0x10000)
        D = 0x10000;
    }
    for (Y=-R*0x10000;Y<=R*0x10000;Y+=D)
    {
      i = (Y+0x8000)/0x10000;
      y = yc+i;
      if ((y<1)||(y>=ys-1))
        continue;
      // left side
      j_out = RR_ii = INTSQRT(R*R-i*i);
      if (xc-j_out<1)
        j_out = xc-1;
      j_inn = (i>0)?i:(-i);
//j_inn = 2*((i>0)?i:(-i));
      if (j_inn<r)
      {
        rr_ii = INTSQRT(r*r-i*i);
        if (j_inn<rr_ii)
          j_inn = rr_ii;
      }
      L = i*i+j_out*j_out;
      currad = INTSQRT(L);
      psim = smoothedim+y*xs+xc-j_out;
      for (j=-j_out;j<=-j_inn;j++)
      {
        // count number of points in circle
        numpoL[currad]++;
        // straight gradient
//        gx = (int)(1.2*(psim[ 1]-psim[-1]));
        gx = 6*(psim[ 1]-psim[-1])/5;
        if (gx<-2)
        {
          gy =psim[ xs]-psim[-xs];
          G = gx*gx+gy*gy;
          if (G>2*36*2)
          {
            // inner product
            S = gx*j+gy*i;
            //
// #ifdef BOTH_GRAD_DIR
            if (((int40)(S*3))*S>((int40)L)*G*2)
              pnProjLN[currad]++;
// #else //BOTH_GRAD_DIR
//            if ((S>0)&&(((int40)(S*3))*S>((int40)L)*G*2))
//             if ((S>0)&&(((int40)(S*4))*S>((int40)L)*G*3))
//if ((S>0)&&(((int40)(S*5))*S>((int40)L)*G*4))
              pnProjLN[currad]++;
// #endif //BOTH_GRAD_DIR
          }
        }
        psim++;
        L -= (-j*2-1);
        if (L<currad*currad)
          currad--;
      }
      // right side
      j_out = RR_ii;
      if (xc+j_out>xs-2)
        j_out = xs-2-xc;
      L = i*i+j_out*j_out;
      currad = INTSQRT(L);
      psim = smoothedim+y*xs+xc+j_out;
      for (j=j_out;j>=j_inn;j--)
      {
        // count number of points in circle
        numpoR[currad]++;
        // straight gradient
//        gx = (int)(1.2*(psim[ 1]-psim[-1]));
        gx = 6*(psim[ 1]-psim[-1])/5;
        if (gx>2)
        {
          gy =psim[ xs]-psim[-xs];
          G = gx*gx+gy*gy;
          if (G>2*36*2)
          {
            // inner product
            S = gx*j+gy*i;
            //
// #ifdef BOTH_GRAD_DIR
            if (((int40)(S*3))*S>((int40)L)*G*2)
              pnProjRN[currad]++;
// #else //BOTH_GRAD_DIR
            //if ((S>0)&&(((int40)(S*3))*S>((int40)L)*G*2))
// if ((S>0)&&(((int40)(S*4))*S>((int40)L)*G*3))
//if ((S>0)&&(((int40)(S*5))*S>((int40)L)*G*4))
              pnProjRN[currad]++;
// #endif //BOTH_GRAD_DIR
          }
        }
        psim--;
        L -= (j*2-1);
        if (L<currad*currad)
          currad--;
      }
    }
    // normalize
    for (i=0;i<r;i++)
    {
      numpoR[i] = numpoR[r];
      numpoL[i] = numpoL[r];
      pnProjRN[i] = pnProjRN[r];
      pnProjLN[i] = pnProjLN[r];
    }
    for (i=0;i<=R;i++)
    {
      pnProjRN[i] = (numpoR[i])?((pnProjRN[i]*1024)/numpoR[i]):0;
      pnProjLN[i] = (numpoL[i])?((pnProjLN[i]*1024)/numpoL[i]):0;
    }
    // process histogram: blur
    IPL_HIST_Blur(numpoR_b,numpoR,R+1,BLUR_HIST_WND);
    IPL_HIST_Blur(numpoL_b,numpoL,R+1,BLUR_HIST_WND);
    IPL_HIST_Blur(pnProjRN_b,pnProjRN,R+1,BLUR_HIST_WND);
    IPL_HIST_Blur(pnProjLN_b,pnProjLN,R+1,BLUR_HIST_WND);
    // correct for small presence
    for (i=1;i<=R;i++)
    {
      L = 10000*numpoR_b[i]/(157*(2*BLUR_HIST_WND+1)*i);
      if (L<20)
        pnProjRN_b[i] = 0;
      L = 10000*numpoL_b[i]/(157*(2*BLUR_HIST_WND+1)*i);
      if (L<20)
        pnProjLN_b[i] = 0;
    }
    pnProjRN_b[0] = pnProjRN_b[1];
    pnProjLN_b[0] = pnProjLN_b[1];
    return ERROR_OK;
}

// same as '5', but with thresholds from '3'
RESULT_CODE IPL_PROJ_FindHoughDonatorProjection6(
  int* pnProjRN_b,          // OUT: circular projection - right side nums
  int* pnProjLN_b,          // OUT: circular projection - left side nums
  const unsigned char* im,  // IN:  image
  int xs,                   // IN:  image size
  int ys,                   // 
  int xc,                   // IN:  ring center
  int yc,                   //
  int r,                    // IN:  inner radius
  int* pR,                  // IN/OUT: outer radius proposed/actually used
  int nMaxRadOfDence,       // IN:  maximum radius for dense processing. 0 - dense always
  void* buf,                // IN:  external buffer
  int* buflen)              // IN/OUT: allocated/used bytes
{
  int i,j,gx,gy,y,G,L,S,*numpoR,*numpoL,currad,R;
  int j_inn,j_out,Y,D,RR_ii,rr_ii;
  const unsigned char* pim;
  int *pnProjRN,*pnProjLN,*numpoR_b,*numpoL_b;

    // check arguments
    if ((buflen==NULL)||(pR==NULL))
      return ERROR_NULL_POINTER;
    // verify radii and calculate max outer radius
    R = *pR;
    if (R<10)
      return ERROR_WRONG_INPUT;
    if ((r+5>R)||(r<5))
      return ERROR_WRONG_INPUT;
    if ((xs<2*r)||(ys<2*r))
      return ERROR_WRONG_INPUT;
    if ((xc-r/2<0)||(xc+r/2>=xs)||(yc-r/2<0)||(yc+r/2>=ys))
      return ERROR_WRONG_INPUT;
    // verify outer radius - bound by maximum span from center to border
    i = xs-1-xc;
    if (i<xc)
      i = xc;
    if (R>i)
      R = i;
    // verify outer radius - bound by image size
    if (2*R>xs)
      R = xs/2;
    // return corrected outer radius
    *pR = R;
    // calc size
    i = (2+4+2)*(R+1)*sizeof(numpoR[0]);
    if (buf==NULL)
    {
      *buflen = i;
      return ERROR_OK;
    }
    if (*buflen<i)
    {
      *buflen = i;
      return ERROR_NULL_POINTER; // ERR_GEN_NO_MEMORY
    }
    // more check arguments
    if ((pnProjRN_b==NULL)||(pnProjLN_b==NULL)||(im==NULL))
      return ERROR_NULL_POINTER;
    /*if ((((ptr_t)im)&3)||(((ptr_t)buf)&3)||(xs&3)||
        (((ptr_t)pnProjRN_b)&3)||(((ptr_t)pnProjLN_b)&3))
      return ERR_GEN_BAD_ALIGNMENT;
    */
    // allocate
    numpoR   = (int*)buf;
    numpoL   = numpoR+(R+1);
    pnProjRN = numpoL+(R+1);
    pnProjLN = pnProjRN+(R+1);
    numpoR_b = pnProjLN+(R+1);
    numpoL_b = numpoR_b+(R+1);
    // clear junk
    memset(pnProjRN,0,sizeof(pnProjRN[0])*(1+R));
    memset(pnProjLN,0,sizeof(pnProjLN[0])*(1+R));
    memset(numpoR,0,sizeof(numpoR[0])*(1+R));
    memset(numpoL,0,sizeof(numpoL[0])*(1+R));
    // along y
    if (!nMaxRadOfDence)
      D = 0x10000;
    else
    {
      D = (r*0x10000)/nMaxRadOfDence;
      if (D<0x10000)
        D = 0x10000;
    }
    for (Y=-R*0x10000;Y<=R*0x10000;Y+=D)
    {
      i = (Y+0x8000)/0x10000;
      y = yc+i;
      if ((y<1)||(y>=ys-1))
        continue;

      // left side
      j_out = RR_ii = INTSQRT(R*R-i*i);
      if (xc-j_out<1)
        j_out = xc-1;
      j_inn = (i>0)?i:(-i);
      if (j_inn<r)
      {
        rr_ii = INTSQRT(r*r-i*i);
        if (j_inn<rr_ii)
          j_inn = rr_ii;
      }
      L = i*i+j_out*j_out;
      currad = INTSQRT(L);
      pim = im+y*xs+xc-j_out;
      for (j=-j_out;j<=-j_inn;j++)
      {
        // count number of points in circle
        numpoL[currad]++;
        // straight gradient
        gx =+pim[xs+1]+2*pim[ 1]+pim[-xs+1]
            -pim[xs-1]-2*pim[-1]-pim[-xs-1];
        if (gx<-4)
        {
          gy =+pim[ xs+1]+2*pim[ xs]+pim[ xs-1]
              -pim[-xs+1]-2*pim[-xs]-pim[-xs-1];
          G = gx*gx+gy*gy;
          if (G>2*36)
//          if (G>36)
          {
            // inner product
            S = gx*j+gy*i;
            //
            if ((S>0)&&(((int40)(S*2))*S>((int40)L)*G))
              pnProjLN[currad]++;
          }
        }
        pim++;
        L -= (-j*2-1);
        if (L<currad*currad)
          currad--;
      }

      // right side
      j_out = RR_ii;
      if (xc+j_out>xs-2)
        j_out = xs-2-xc;
      L = i*i+j_out*j_out;
      currad = INTSQRT(L);
      pim = im+y*xs+xc+j_out;
      for (j=j_out;j>=j_inn;j--)
      {
        // count number of points in circle
        numpoR[currad]++;
        // straight gradient
        gx =+pim[xs+1]+2*pim[ 1]+pim[-xs+1]
            -pim[xs-1]-2*pim[-1]-pim[-xs-1];
        if (gx>4)
        {
          gy =+pim[ xs+1]+2*pim[ xs]+pim[ xs-1]
              -pim[-xs+1]-2*pim[-xs]-pim[-xs-1];
          G = gx*gx+gy*gy;
          if (G>2*36)
//          if (G>36)
          {
            // inner product
            S = gx*j+gy*i;
            //
            if ((S>0)&&(((int40)(S*2))*S>((int40)L)*G))
              pnProjRN[currad]++;
          }
        }
        pim--;
        L -= (j*2-1);
        if (L<currad*currad)
          currad--;
      }
    }
    // normalize
    for (i=0;i<r;i++)
    {
      numpoR[i] = numpoR[r];
      numpoL[i] = numpoL[r];
      pnProjRN[i] = pnProjRN[r];
      pnProjLN[i] = pnProjLN[r];
    }
    for (i=0;i<=R;i++)
    {
      pnProjRN[i] = (numpoR[i])?((pnProjRN[i]*1024)/numpoR[i]):0;
      pnProjLN[i] = (numpoL[i])?((pnProjLN[i]*1024)/numpoL[i]):0;
    }
    // process histogram: blur
    IPL_HIST_Blur(numpoR_b,numpoR,R+1,BLUR_HIST_WND);
    IPL_HIST_Blur(numpoL_b,numpoL,R+1,BLUR_HIST_WND);
    IPL_HIST_Blur(pnProjRN_b,pnProjRN,R+1,BLUR_HIST_WND);
    IPL_HIST_Blur(pnProjLN_b,pnProjLN,R+1,BLUR_HIST_WND);
    // correct for small presence
    for (i=1;i<=R;i++)
    {
      L = 10000*numpoR_b[i]/(157*(2*BLUR_HIST_WND+1)*i);
      if (L<20)
        pnProjRN_b[i] = 0;
      L = 10000*numpoL_b[i]/(157*(2*BLUR_HIST_WND+1)*i);
      if (L<20)
        pnProjLN_b[i] = 0;
    }
    pnProjRN_b[0] = pnProjRN_b[1];
    pnProjLN_b[0] = pnProjLN_b[1];
    return ERROR_OK;
}

// calculate left and right side projection histograms in a concentric ring
RESULT_CODE IPL_PROJ_FindHoughDonatorProjection5(
  int* pnProjRN_b,          // OUT: circular projection - right side nums
  int* pnProjLN_b,          // OUT: circular projection - left side nums
  const unsigned char* im,  // IN:  image
  int xs,                   // IN:  image size
  int ys,                   // 
  int xc,                   // IN:  ring center
  int yc,                   //
  int r,                    // IN:  inner radius
  int* pR,                  // IN/OUT: outer radius proposed/actually used
  int nMaxRadOfDence,       // IN:  maximum radius for dense processing. 0 - dense always
  void* buf,                // IN:  external buffer
  int* buflen)              // IN/OUT: allocated/used bytes
{
  int i,j,gx,gy,y,G,L,S,*numpoR,*numpoL,currad,R;
  int j_inn,j_out,Y,D,RR_ii,rr_ii;
  const unsigned char* pim;
  int *pnProjRN,*pnProjLN,*numpoR_b,*numpoL_b;

    // check arguments
    if ((buflen==NULL)||(pR==NULL))
      return ERROR_NULL_POINTER;
    // verify radii and calculate max outer radius
    R = *pR;
    if (R<10)
      return ERROR_WRONG_INPUT;
    if ((r+5>R)||(r<5))
      return ERROR_WRONG_INPUT;
    if ((xs<2*r)||(ys<2*r))
      return ERROR_WRONG_INPUT;
    if ((xc-r/2<0)||(xc+r/2>=xs)||(yc-r/2<0)||(yc+r/2>=ys))
      return ERROR_WRONG_INPUT;
    // verify outer radius - bound by maximum span from center to border
    i = xs-1-xc;
    if (i<xc)
      i = xc;
    if (R>i)
      R = i;
    // verify outer radius - bound by image size
    if (2*R>xs)
      R = xs/2;
    // return corrected outer radius
    *pR = R;
    // calc size
    i = (2+4+2)*(R+1)*sizeof(numpoR[0]);
    if (buf==NULL)
    {
      *buflen = i;
      return ERROR_OK;
    }
    if (*buflen<i)
    {
      *buflen = i;
      return ERROR_MEMORY;
    }
    // more check arguments
    if ((pnProjRN_b==NULL)||(pnProjLN_b==NULL)||(im==NULL))
      return ERROR_NULL_POINTER;
    
    /*
    if ((((ptr_t)im)&3)||(((ptr_t)buf)&3)||(xs&3)||
        (((ptr_t)pnProjRN_b)&3)||(((ptr_t)pnProjLN_b)&3))
      return ERR_GEN_BAD_ALIGNMENT;
    */

    // allocate
    numpoR   = (int*)buf;
    numpoL   = numpoR+(R+1);
    pnProjRN = numpoL+(R+1);
    pnProjLN = pnProjRN+(R+1);
    numpoR_b = pnProjLN+(R+1);
    numpoL_b = numpoR_b+(R+1);
    // clear junk
    memset(pnProjRN,0,sizeof(pnProjRN[0])*(1+R));
    memset(pnProjLN,0,sizeof(pnProjLN[0])*(1+R));
    memset(numpoR,0,sizeof(numpoR[0])*(1+R));
    memset(numpoL,0,sizeof(numpoL[0])*(1+R));
    // along y
    if (!nMaxRadOfDence)
      D = 0x10000;
    else
    {
      D = (r*0x10000)/nMaxRadOfDence;
      if (D<0x10000)
        D = 0x10000;
    }
    for (Y=-R*0x10000;Y<=R*0x10000;Y+=D)
    {
      i = (Y+0x8000)/0x10000;
      y = yc+i;
      if ((y<1)||(y>=ys-1))
        continue;
      // left side
      j_out = RR_ii = INTSQRT(R*R-i*i);
      if (xc-j_out<1)
        j_out = xc-1;
      j_inn = (i>0)?i:(-i);
      if (j_inn<r)
      {
        rr_ii = INTSQRT(r*r-i*i);
        if (j_inn<rr_ii)
          j_inn = rr_ii;
      }
      L = i*i+j_out*j_out;
      currad = INTSQRT(L);
      pim = im+y*xs+xc-j_out;
      for (j=-j_out;j<=-j_inn;j++)
      {
        // count number of points in circle
        numpoL[currad]++;
        // straight gradient
        gx =+pim[xs+1]+2*pim[ 1]+pim[-xs+1]
            -pim[xs-1]-2*pim[-1]-pim[-xs-1];
        if (gx<-6)
        {
          gy =+pim[ xs+1]+2*pim[ xs]+pim[ xs-1]
              -pim[-xs+1]-2*pim[-xs]-pim[-xs-1];
          G = gx*gx+gy*gy;
          if (G>2*64)
          {
            // inner product
            S = gx*j+gy*i;
            //
            if ((S>0)&&(((int40)(S*2))*S>((int40)L)*G))
              pnProjLN[currad]++;
          }
        }
        pim++;
        L -= (-j*2-1);
        if (L<currad*currad)
          currad--;
      }
      // right side
      j_out = RR_ii;
      if (xc+j_out>xs-2)
        j_out = xs-2-xc;
      L = i*i+j_out*j_out;
      currad = INTSQRT(L);
      pim = im+y*xs+xc+j_out;
      for (j=j_out;j>=j_inn;j--)
      {
        // count number of points in circle
        numpoR[currad]++;
        // straight gradient
        gx =+pim[xs+1]+2*pim[ 1]+pim[-xs+1]
            -pim[xs-1]-2*pim[-1]-pim[-xs-1];
        if (gx>6)
        {
          gy =+pim[ xs+1]+2*pim[ xs]+pim[ xs-1]
              -pim[-xs+1]-2*pim[-xs]-pim[-xs-1];
          G = gx*gx+gy*gy;
          if (G>2*64)
          {
            // inner product
            S = gx*j+gy*i;
            //
            if ((S>0)&&(((int40)(S*2))*S>((int40)L)*G))
              pnProjRN[currad]++;
          }
        }
        pim--;
        L -= (j*2-1);
        if (L<currad*currad)
          currad--;
      }
    }
    // normalize
    for (i=0;i<r;i++)
    {
      numpoR[i] = numpoR[r];
      numpoL[i] = numpoL[r];
      pnProjRN[i] = pnProjRN[r];
      pnProjLN[i] = pnProjLN[r];
    }
    for (i=0;i<=R;i++)
    {
      pnProjRN[i] = (numpoR[i])?((pnProjRN[i]*1024)/numpoR[i]):0;
      pnProjLN[i] = (numpoL[i])?((pnProjLN[i]*1024)/numpoL[i]):0;
    }
    // process histogram: blur
    IPL_HIST_Blur(numpoR_b,numpoR,R+1,BLUR_HIST_WND);
    IPL_HIST_Blur(numpoL_b,numpoL,R+1,BLUR_HIST_WND);
    IPL_HIST_Blur(pnProjRN_b,pnProjRN,R+1,BLUR_HIST_WND);
    IPL_HIST_Blur(pnProjLN_b,pnProjLN,R+1,BLUR_HIST_WND);
    // correct for small presence
    for (i=1;i<=R;i++)
    {
      L = 10000*numpoR_b[i]/(157*(2*BLUR_HIST_WND+1)*i);
      if (L<20)
        pnProjRN_b[i] = 0;
      L = 10000*numpoL_b[i]/(157*(2*BLUR_HIST_WND+1)*i);
      if (L<20)
        pnProjLN_b[i] = 0;
    }
    pnProjRN_b[0] = pnProjRN_b[1];
    pnProjLN_b[0] = pnProjLN_b[1];
    return ERROR_OK;
}

// calculate left and right side projection histograms in a concentric ring
RESULT_CODE IPL_PROJ_FindHoughDonatorProjection5_Mask(
  int* pnProjRN_b,          // OUT: circular projection - right side nums
  int* pnProjLN_b,          // OUT: circular projection - left side nums
  const unsigned char* im,  // IN:  image
  const unsigned char* mask,// IN:  image
  int xs,                   // IN:  image size
  int ys,                   // 
  int xc,                   // IN:  ring center
  int yc,                   //
  int r,                    // IN:  inner radius
  int* pR,                  // IN/OUT: outer radius proposed/actually used
  int nMaxRadOfDence,       // IN:  maximum radius for dense processing. 0 - dense always
  void* buf,                // IN:  external buffer
  int* buflen)              // IN/OUT: allocated/used bytes
{
  int i,j,gx,gy,y,G,L,S,*numpoR,*numpoL,currad,R;
  int j_inn,j_out,Y,D,RR_ii,rr_ii;
  const unsigned char* pim;
  int *pnProjRN,*pnProjLN,*numpoR_b,*numpoL_b;

    // check arguments
    if ((buflen==NULL)||(pR==NULL))
      return ERROR_NULL_POINTER;
    // verify radii and calculate max outer radius
    R = *pR;
    if (R<10)
      return ERROR_WRONG_INPUT;
    if ((r+5>R)||(r<5))
      return ERROR_WRONG_INPUT;
    if ((xs<2*r)||(ys<2*r))
      return ERROR_WRONG_INPUT;
    if ((xc-r/2<0)||(xc+r/2>=xs)||(yc-r/2<0)||(yc+r/2>=ys))
      return ERROR_WRONG_INPUT;
    // verify outer radius - bound by maximum span from center to border
    i = xs-1-xc;
    if (i<xc)
      i = xc;
    if (R>i)
      R = i;
    // verify outer radius - bound by image size
    if (2*R>xs)
      R = xs/2;
    // return corrected outer radius
    *pR = R;
    // calc size
    i = (2+4+2)*(R+1)*sizeof(numpoR[0]);
    if (buf==NULL)
    {
      *buflen = i;
      return ERROR_OK;
    }
    if (*buflen<i)
    {
      *buflen = i;
      return ERROR_MEMORY;
    }
    // more check arguments
    if ((pnProjRN_b==NULL)||(pnProjLN_b==NULL)||(im==NULL))
      return ERROR_NULL_POINTER;
    /*if ((((ptr_t)im) & 3) || (((ptr_t)buf) & 3) || (xs & 3) ||
        (((ptr_t)pnProjRN_b) & 3) || (((ptr_t)pnProjLN_b) & 3))
        return ERROR_MEMORY; // ERR_GEN_BAD_ALIGNMENT;
    */
    // allocate
    numpoR   = (int*)buf;
    numpoL   = numpoR+(R+1);
    pnProjRN = numpoL+(R+1);
    pnProjLN = pnProjRN+(R+1);
    numpoR_b = pnProjLN+(R+1);
    numpoL_b = numpoR_b+(R+1);
    // clear junk
    memset(pnProjRN,0,sizeof(pnProjRN[0])*(1+R));
    memset(pnProjLN,0,sizeof(pnProjLN[0])*(1+R));
    memset(numpoR,0,sizeof(numpoR[0])*(1+R));
    memset(numpoL,0,sizeof(numpoL[0])*(1+R));
    // along y
    if (!nMaxRadOfDence)
      D = 0x10000;
    else
    {
      D = (r*0x10000)/nMaxRadOfDence;
      if (D<0x10000)
        D = 0x10000;
    }
    for (Y=-R*0x10000;Y<=R*0x10000;Y+=D)
    {
      i = (Y+0x8000)/0x10000;
      y = yc+i;
      if ((y<1)||(y>=ys-1))
        continue;
      // left side
      j_out = RR_ii = INTSQRT(R*R-i*i);
      if (xc-j_out<1)
        j_out = xc-1;
      j_inn = (i>0)?i:(-i);
      if (j_inn<r)
      {
        rr_ii = INTSQRT(r*r-i*i);
        if (j_inn<rr_ii)
          j_inn = rr_ii;
      }
      L = i*i+j_out*j_out;
      currad = INTSQRT(L);
      pim = im+y*xs+xc-j_out;
      for (j=-j_out;j<=-j_inn;j++)
      {
        if (mask[pim-im]==0)
        {
          // count number of points in circle
          numpoL[currad]++;
          // straight gradient
          gx =+pim[xs+1]+2*pim[ 1]+pim[-xs+1]
              -pim[xs-1]-2*pim[-1]-pim[-xs-1];
          if (gx<-6)
          {
            gy =+pim[ xs+1]+2*pim[ xs]+pim[ xs-1]
                -pim[-xs+1]-2*pim[-xs]-pim[-xs-1];
            G = gx*gx+gy*gy;
            if (G>2*64)
            {
              // inner product
              S = gx*j+gy*i;
              //
              if ((S>0)&&(((int40)(S*2))*S>((int40)L)*G))
                pnProjLN[currad]++;
            }
          }
        }
        pim++;
        L -= (-j*2-1);
        if (L<currad*currad)
          currad--;
      }
      // right side
      j_out = RR_ii;
      if (xc+j_out>xs-2)
        j_out = xs-2-xc;
      L = i*i+j_out*j_out;
      currad = INTSQRT(L);
      pim = im+y*xs+xc+j_out;
      for (j=j_out;j>=j_inn;j--)
      {
        if (mask[pim-im]==0)
        {
          // count number of points in circle
          numpoR[currad]++;
          // straight gradient
          gx =+pim[xs+1]+2*pim[ 1]+pim[-xs+1]
              -pim[xs-1]-2*pim[-1]-pim[-xs-1];
          if (gx>6)
          {
            gy =+pim[ xs+1]+2*pim[ xs]+pim[ xs-1]
                -pim[-xs+1]-2*pim[-xs]-pim[-xs-1];
            G = gx*gx+gy*gy;
            if (G>2*64)
            {
              // inner product
              S = gx*j+gy*i;
              //
              if ((S>0)&&(((int40)(S*2))*S>((int40)L)*G))
                pnProjRN[currad]++;
            }
          }
        }
        pim--;
        L -= (j*2-1);
        if (L<currad*currad)
          currad--;
      }
    }
    // normalize
    for (i=0;i<r;i++)
    {
      numpoR[i] = numpoR[r];
      numpoL[i] = numpoL[r];
      pnProjRN[i] = pnProjRN[r];
      pnProjLN[i] = pnProjLN[r];
    }
    for (i=0;i<=R;i++)
    {
      pnProjRN[i] = (numpoR[i])?((pnProjRN[i]*1024)/numpoR[i]):0;
      pnProjLN[i] = (numpoL[i])?((pnProjLN[i]*1024)/numpoL[i]):0;
    }
    // process histogram: blur
    IPL_HIST_Blur(numpoR_b,numpoR,R+1,BLUR_HIST_WND);
    IPL_HIST_Blur(numpoL_b,numpoL,R+1,BLUR_HIST_WND);
    IPL_HIST_Blur(pnProjRN_b,pnProjRN,R+1,BLUR_HIST_WND);
    IPL_HIST_Blur(pnProjLN_b,pnProjLN,R+1,BLUR_HIST_WND);
    // correct for small presence
    for (i=1;i<=R;i++)
    {
      L = 10000*numpoR_b[i]/(157*(2*BLUR_HIST_WND+1)*i);
      if (L<20)
        pnProjRN_b[i] = 0;
      L = 10000*numpoL_b[i]/(157*(2*BLUR_HIST_WND+1)*i);
      if (L<20)
        pnProjLN_b[i] = 0;
    }
    pnProjRN_b[0] = pnProjRN_b[1];
    pnProjLN_b[0] = pnProjLN_b[1];
    return ERROR_OK;
}

static RESULT_CODE ownIVIR_SelectPupir(
    char* casechar,   // OUT: character for diagnose
    int* pXpl,
    int* pXpr,
    int* pQp,
    int* pXil,
    int* pXir,
    int* pQi,
    const int* projR,
    const int* anMaxsR,
    int nMaxR,
    const int* projL,
    const int* anMaxsL,
    int nMaxL,
    int projlen,
    const unsigned char* im,
    int xs,
    int ys,
    int xc,
    int yc,
    int mode,
    int br_size,
    int br_spread)
{
    int Xpl, Xpr, Xil, Xir, i, j, k, l, rp, ri, Nr, Nl;
    int Qs, Qb, Qd, Qv, Q, max3r, max3l, zs, maxXpl, maxXpr, maxXil, maxXir;
    int valXpl, valXpr, valXil, valXir, PosMaxR, PosMaxL, isLeft = 0, isRight = 0;
    SBazRadHyp* psbrh = NULL;
    int nbrh_used = 0, nbrh_alloced = 0;

    // set side information
    if (*casechar == 'L')
        isLeft = 1;
    if (*casechar == 'R')
        isRight = 1;
    zs = (xs < ys) ? xs : ys;
    // default - not found
    *pXpl = -1;
    *pXpr = -1;
    *pXil = -1;
    *pXir = -1;
    *casechar = '-';
    *pQp = *pQi = -1;
    // init positions
    maxXpl = maxXpr = maxXil = maxXir = -1;
    if ((nMaxR < 1) || (nMaxL < 1))
    { // no peaks found at least at one side (virtually impossible)
        *casechar = 'N';
        goto quita;
    }
    // find position of max in left and right
    PosMaxR = anMaxsR[0];
    for (i = 1; i < nMaxR; i++)
        if (projR[PosMaxR] < projR[anMaxsR[2 * i]])
            PosMaxR = anMaxsR[2 * i];
    PosMaxL = anMaxsL[0];
    for (i = 1; i < nMaxL; i++)
        if (projL[PosMaxL] < projL[anMaxsL[2 * i]])
            PosMaxL = anMaxsL[2 * i];
    if ((nMaxR < 2) || (nMaxL < 2))
        // only one peak found (very low probability of this)
        goto quita;
    // calculate sum of maxima
    //    for (Qs=i=0;i<nMaxR;i++)
    //    Qs += projR[anMaxsR[i*2]];
    //Qs /= nMaxR;
    for (Qv = i = 0; i < nMaxR - 1; i++)
        Qv += projR[anMaxsR[i * 2 + 1]];
    Qv /= (nMaxR - 1);
    if (2 * Qv > projR[PosMaxR])
    {
        *casechar = 'X';
        goto quita;
    }
    
    for (Qv = i = 0; i < nMaxL - 1; i++)
        Qv += projL[anMaxsL[i * 2 + 1]];
    Qv /= (nMaxL - 1);
    if (2 * Qv > projL[PosMaxL])
    {
        *casechar = 'X';
        goto quita;
    }
    // two or more peaks in each side
    for (i = 0; i < nMaxR - 1; i++)
    {
        Xpr = anMaxsR[2 * i];
        for (j = i + 1; j < nMaxR; j++)
        {
            Xir = anMaxsR[2 * j];
            for (k = 0; k < nMaxL - 1; k++)
            {
                Xpl = anMaxsL[2 * k];
                for (l = k + 1; l < nMaxL; l++)
                {
                    Xil = anMaxsL[2 * l];
                    // use side information
                    if (isRight && (!(Xir - Xpr < Xil - Xpl)))
                        continue;
                    if (isLeft && (!(Xir - Xpr > Xil - Xpl)))
                        continue;
                    // use predefinitions
                    if (br_size > 0)
                        if ((Xil + Xir < 2 * (br_size - br_spread)) || (Xil + Xir > 2 * (br_size + br_spread)))
                            continue;
                    // at least one of maxima should be sufficiently big
                    if (((projL[PosMaxL] > 2 * projL[Xpl]) && (projL[PosMaxL] > 4 * projL[Xil])) ||
                        ((projR[PosMaxR] > 2 * projR[Xpr]) && (projR[PosMaxR] > 4 * projR[Xir])))
                        continue;
                    rp = (Xpl + Xpr) / 2;
                    ri = (Xil + Xir) / 2;
                    // filters
                    if (rp * 7 < ri)  // mia120427
                        continue;   // pupil should be at least 1/7 of iris
                    if (rp * 4 > ri * 3)
                        continue;   // pupil should be at most 3/4 of iris
                    if ((2 * (Xil - Xpl) < (Xir - Xpr)) || (2 * (Xir - Xpr) < (Xil - Xpl)))
                        continue;   // iris sides should not differ strongly
                    if ((ri < zs / 10) || (ri > zs / 2) || (rp < 5) || (rp > zs / 3))
                        continue;   // reject too big or too small eyes relative to image size
                                    // calculate peak height ratios
                    valXpl = projL[Xpl];
                    valXpr = projR[Xpr];
                    valXil = projL[Xil];
                    valXir = projR[Xir];
                    //== qualities
                    // - pupil/iris decentration
                    if (Xir - Xpr < Xil - Xpl)
                        Qd = 100 * (Xir - Xpr) / (Xil - Xpl);
                    else
                        Qd = 100 * (Xil - Xpl) / (Xir - Xpr);
                    if (Qd < 55)
                        Qd = 0;
                    else
                        if (Qd < 67)
                            Qd = (Qd - 55)*(100 - 0) / (67 - 55) + 0;
                        else
                            Qd = 100;
                    // - center/iris decentration
                    Qb = 50 * (Xir - Xil) / (Xil + Xir);
                    if (Qb < 0)
                        Qb = -Qb;
                    // - pupil/iris size
                    Qs = 100 * (Xpl + Xpr) / (Xil + Xir);
                    if (Qs < 20)
                        Qs = (Qs - 14)*(100 - 50) / (20 - 14) + 50;
                    else
                        if (Qs < 70)
                            Qs = 100;
                        else
                            Qs = (Qs - 70)*(50 - 100) / (75 - 70) + 100;
                    
                    Qb = 100 - Qb;
                    // -value
                    Qv = 50 * (valXpr + valXpl + valXir + valXil) / (projL[PosMaxL] + projR[PosMaxR]);
                    //            Q = Qv*Qd*Qb/100+1;
                    Q = Qv*Qd*Qb*Qs / 10000 + 1;
                    if (nbrh_used == nbrh_alloced)
                        psbrh = (SBazRadHyp*)realloc(psbrh, (nbrh_alloced += 1024) * sizeof(psbrh[0]));
                    psbrh[nbrh_used].q = Q;
                    psbrh[nbrh_used].xpl = Xpl;
                    psbrh[nbrh_used].xpr = Xpr;
                    psbrh[nbrh_used].xil = Xil;
                    psbrh[nbrh_used].xir = Xir;
                    psbrh[nbrh_used].vpl = valXpl;
                    psbrh[nbrh_used].vpr = valXpr;
                    psbrh[nbrh_used].vil = valXil;
                    psbrh[nbrh_used].vir = valXir;
                    psbrh[nbrh_used].casechar = '0';
                    nbrh_used++;
                }
            }
        }
    }
    //bongo:;
    // hypotheses are enumerated. Now choose. 
    if (!nbrh_used)
        // no hypothesis at all - nothing found
        goto quita;
    if (nbrh_used == 1)
    { // only one single hypothesis 
        *casechar = 's';//psbrh[0].casechar+4;
        *pXpl = psbrh[0].xpl;
        *pXpr = psbrh[0].xpr;
        *pXil = psbrh[0].xil;
        *pXir = psbrh[0].xir;
        *pQp = *pQi = 50;
        goto quita;
    }
    qsort(psbrh, nbrh_used, sizeof(psbrh[0]), ownIVIR_sortSBRH);
    // calc ratio of pupil and iris peak heights
    Nr = psbrh[0].vpr * 1024 / (psbrh[0].vir + 1);
    Nl = psbrh[0].vpl * 1024 / (psbrh[0].vil + 1);
    if ((Nr > 1024 / 3) && (Nr < 1024 * 3) && (Nl > 1024 / 3) && (Nl < 1024 * 3) &&
        (1000 * psbrh[1].q <= 900 * psbrh[0].q))
    { // almost guarantied
        *casechar = psbrh[0].casechar;
        *pXpl = psbrh[0].xpl;
        *pXpr = psbrh[0].xpr;
        *pXil = psbrh[0].xil;
        *pXir = psbrh[0].xir;
        Q = 100;
        *pQp = *pQi = Q;
        goto quita;
    }
    if ((Nr > 800) && (Nr < 1311) && (Nl > 800) && (Nl < 1311))
    { // almost guarantied - very good pupil and iris
        *casechar = psbrh[0].casechar + 1;
        *pXpl = psbrh[0].xpl;
        *pXpr = psbrh[0].xpr;
        *pXil = psbrh[0].xil;
        *pXir = psbrh[0].xir;
        Q = 100;
        *pQp = *pQi = Q;
        goto quita;
    }
    if ((Nr > 1024 / 8) && (Nr < 1024 * 8) && (Nl > 1024 / 8) && (Nl < 1024 * 8) &&
        (1000 * psbrh[1].q <= 915 * psbrh[0].q))
    { // bit worse
        *casechar = psbrh[0].casechar + 2;
        *pXpl = psbrh[0].xpl;
        *pXpr = psbrh[0].xpr;
        *pXil = psbrh[0].xil;
        *pXir = psbrh[0].xir;
        Q = 95;
        *pQp = *pQi = Q;
        goto quita;
    }
    if ((Nr > 1024 / 2) && (Nr < 1024 * 2) && (Nl > 1024 / 2) && (Nl < 1024 * 2) &&
        (1000 * psbrh[1].q <= 930 * psbrh[0].q))
    { // worse
        *casechar = psbrh[0].casechar + 3;
        *pXpl = psbrh[0].xpl;
        *pXpr = psbrh[0].xpr;
        *pXil = psbrh[0].xil;
        *pXir = psbrh[0].xir;
        Q = 90;
        *pQp = *pQi = Q;
        goto quita;
    }
    //    if ((Nr>1024/4)&&(Nr<1024*4)&&(Nl>1024/4)&&(Nl<1024*4)&&
    //      (1000*psbrh[1].q<=930*psbrh[0].q))
    if ((Nr > 1024 / 2) && (Nr < 1024 * 2) && (Nl > 1024 / 2) && (Nl < 1024 * 2))
    { // more worse
        *casechar = psbrh[0].casechar + 4;
        *pXpl = psbrh[0].xpl;
        *pXpr = psbrh[0].xpr;
        *pXil = psbrh[0].xil;
        *pXir = psbrh[0].xir;
        Q = 85;
        *pQp = *pQi = Q;
        goto quita;
    }
    if ((Nr > 1024 / 3) && (Nr < 1024 * 3) && (Nl > 1024 / 3) && (Nl < 1024 * 3) &&
        ((psbrh[1].xpl == psbrh[0].xpl) && (psbrh[1].xpr == psbrh[0].xpr) &&
        ((psbrh[1].xil == psbrh[0].xil) || (psbrh[1].xir == psbrh[0].xir))))
    { // max and submax pupils match
        *casechar = 'A';
        *pXpl = psbrh[0].xpl;
        *pXpr = psbrh[0].xpr;
        Q = 75;
        *pQp = Q;
        *pQi = 0;
        goto quita;
    }
    if ((Nr > 1024 / 7) && (Nr < 1024 * 7) && (Nl > 1024 / 7) && (Nl < 1024 * 7) &&
        ((psbrh[1].xil == psbrh[0].xil) && (psbrh[1].xir == psbrh[0].xir) &&
        ((psbrh[1].xpl == psbrh[0].xpl) || (psbrh[1].xpr == psbrh[0].xpr))))
    { // max and submax irises match
        *casechar = 'B';
        *pXil = psbrh[0].xil;
        *pXir = psbrh[0].xir;
        Q = 80;
        *pQi = Q;
        *pQp = 0;
        goto quita;
    }
    // just make a silly guess
    ownIVIR_GetMaxSubMax(&maxXpl, &maxXil, projL, anMaxsL, nMaxL, &max3l);
    ownIVIR_GetMaxSubMax(&maxXpr, &maxXir, projR, anMaxsR, nMaxR, &max3r);
    if (maxXpl < 0)
        maxXpl = 0;
    if (maxXpr < 0)
        maxXpr = 0;
    if ((ownIVIR_difil(anMaxsL[2 * maxXil], anMaxsR[2 * maxXir]) < 150) &&
        (3 * projL[anMaxsL[2 * maxXil]] > 2 * projL[PosMaxL]) &&
        (3 * projR[anMaxsR[2 * maxXir]] > 2 * projR[PosMaxR]))
    { // second maxima positions collimate, both maxima values big => iris
        if (anMaxsL[2 * maxXil] + anMaxsR[2 * maxXir] > 75 * 2)
        {
            *casechar = 'I';
            *pXil = anMaxsL[2 * maxXil];
            *pXir = anMaxsR[2 * maxXir];
            Q = (mode&IVIR_ARRPOXIR_PRPRPR) ? 45 : 83;
            *pQi = Q;
            *pQp = 0;
        }
        else
        {
            *casechar = 'J';
            *pXpl = anMaxsL[2 * maxXil];
            *pXpr = anMaxsR[2 * maxXir];
            //      Q = 100-400*submaxQ/maxQ;
            //    if (Q<=0)
            Q = (mode&IVIR_ARRPOXIR_PRPRPR) ? 10 : 22;
            *pQi = 0;
            *pQp = Q;
        }
        goto quita;
    }
    if ((ownIVIR_difil(anMaxsL[2 * maxXpl], anMaxsR[2 * maxXpr]) < 150) &&
        (4 * projL[anMaxsL[2 * maxXpl]] >= 2 * projL[PosMaxL]) &&
        (4 * projR[anMaxsR[2 * maxXpr]] >= 2 * projR[PosMaxR]))
    { // first maxima positions collimate, both maxima values big =>
      // this is one circle, but pupil or iris?
        int inL, outL, inR, outR, HL, HR;

        for (inL = i = 0; i < maxXpl; i++)
            inL += projL[anMaxsL[i * 2]];
        for (inR = i = 0; i < maxXpr; i++)
            inR += projR[anMaxsR[i * 2]];
        for (outL = 0, i = maxXpl + 1; i < nMaxL; i++)
            outL += projL[anMaxsL[i * 2]];
        for (outR = 0, i = maxXpr + 1; i < nMaxR; i++)
            outR += projR[anMaxsR[i * 2]];
        HL = projL[anMaxsL[maxXpl * 2]];
        HR = projR[anMaxsR[maxXpr * 2]];
        if ((4 * (outL + outR) <= 5 * (HL + HR)) && (inR + inL >= (HL + HR) / 5))
            //      if (0)
        {
            *casechar = 'X';
            *pXil = anMaxsL[2 * maxXpl];
            *pXir = anMaxsR[2 * maxXpr];
            Q = 45;
            *pQi = Q;
            *pQp = 0;
        }
        else
        { // sum of maxima inside is too small, or outsize too big => pupil
            *pXpl = anMaxsL[2 * maxXpl] - 2;
            *pXpr = anMaxsR[2 * maxXpr] - 2;
            if (4 * (outL + outR) <= 5 * (HL + HR))
            {
                *casechar = 'Z';
                Q = (mode&IVIR_ARRPOXIR_PRPRPR) ? 55 : 70;
            }
            else
            {
                *casechar = 'z';
                Q = (mode&IVIR_ARRPOXIR_PRPRPR) ? 50 : 25;
            }
            *pQp = Q;
            *pQi = 0;
        }
        goto quita;
    }
    if ((ownIVIR_difil(anMaxsL[2 * maxXpl], anMaxsR[2 * maxXir]) < 100) &&
        (3 * projL[anMaxsL[2 * maxXpl]] > 2 * projL[PosMaxL]) &&
        (3 * projR[anMaxsR[2 * maxXir]] > 2 * projR[PosMaxR]))
    { // consider as iris
        j = anMaxsL[2 * maxXpl];
        i = anMaxsR[2 * maxXir];
        // !!! size hack - decision depends on image size !!!
        if (((i < projlen / 4) && (j < projlen / 4)) || (i + j < 75))
        {
            Q = (mode&IVIR_ARRPOXIR_PRPRPR) ? 35 : 20;
            *casechar = 'U';
            *pXpr = i - 2;
            *pXpl = j - 2;
            *pQp = Q;
            *pQi = 0;
        }
        else
        {
            Q = (mode&IVIR_ARRPOXIR_PRPRPR) ? 30 : 20;
            *casechar = 'u';
            *pXir = i;
            *pXil = j;
            *pQi = Q;
            *pQp = 0;
        }
        goto quita;
    }
    if ((ownIVIR_difil(anMaxsL[2 * maxXil], anMaxsR[2 * maxXpr]) < 100) &&
        (3 * projL[anMaxsL[2 * maxXil]] > 2 * projL[PosMaxL]) &&
        (3 * projR[anMaxsR[2 * maxXpr]] > 2 * projR[PosMaxR]))
    { // consider as iris
        j = anMaxsL[2 * maxXil];
        i = anMaxsR[2 * maxXpr];
        // !!! size hack - decision depends on image size !!!
        if ((i < projlen / 4) && (j < projlen / 4))
        {
            Q = (mode&IVIR_ARRPOXIR_PRPRPR) ? 40 : 10;
            *casechar = 'V';
            *pXpr = i - 2;
            *pXpl = j - 2;
            *pQp = Q;
            *pQi = 0;
        }
        else
        {
            Q = (mode&IVIR_ARRPOXIR_PRPRPR) ? 25 : 10;
            *casechar = 'v';
            *pXir = i;
            *pXil = j;
            *pQi = Q;
            *pQp = 0;
        }
        goto quita;
    }
    i = anMaxsR[2 * ((projR[anMaxsR[2 * maxXpr]] > projR[anMaxsR[2 * maxXir]]) ? maxXpr : maxXir)];
    j = anMaxsL[2 * ((projL[anMaxsL[2 * maxXpl]] > projL[anMaxsL[2 * maxXil]]) ? maxXpl : maxXil)];
    // !!! size hack - decision depends on image size !!!
    if ((i < 2 * 3 * projlen / (5 * 4)) && (j < 2 * 3 * projlen / (5 * 4)))
    {
        *casechar = 'g';
        *pXpr = i - 2;
        *pXpl = j - 2;
        *pQp = (mode&IVIR_ARRPOXIR_PRPRPR) ? 20 : 15;
        *pQi = 0;
    }
    else
    {
        *casechar = 'G';
        *pXir = i;
        *pXil = j;
        *pQi = (mode&IVIR_ARRPOXIR_PRPRPR) ? 15 : 2;
        *pQp = 0;
    }
quita:
    if ((*casechar >= '0') && (*casechar <= '4'))
    { // try a correction
        i = ownIVIR_CorrectBazradHyp(0, psbrh, nbrh_used,
            *pXpl, *pXpr, *pXil, *pXir, projR, anMaxsR, nMaxR, projL, anMaxsL, nMaxL, projlen);
        j = ownIVIR_CorrectBazradHyp(1, psbrh, nbrh_used,
            *pXpl, *pXpr, *pXil, *pXir, projR, anMaxsR, nMaxR, projL, anMaxsL, nMaxL, projlen);
        if ((i > 0) && (j > 0))
            *casechar = '#';
        if ((i > 0) && (j <= 0))
        {
            *pXir = i;
            *casechar = '@';
        }
        if ((i <= 0) && (j > 0))
        {
            *pXil = j;
            *casechar = '@';
        }
    }
    if (psbrh)
        free(psbrh);
    return ERROR_OK;
}


// detect approximate pupil and iris by projection method
RESULT_CODE IVIR_PrPu(
    SPupilInfo* psPI,         // OUT: pupil data
    SIrisInfo* psII,          // OUT: iris data
    const SCenterInfo* psCI,  // IN:  center data
    const unsigned char* im,  // IN:  source image
    int xs,                   // IN:  image size
    int ys,                   //
    int mode, // mask: 1- draw projs to buf, 2- use projs from buf
    int br_size,    // predefined iris size (-1 - not predefined)
    int br_spread,  // predefined spread of iris size (-1 - not predefined)
    int* projR,
    int* projL,
    const char* nambeg)       // used for debug, set to NULL
{
    RESULT_CODE res = ERROR_OK;
    int sz, r, R, szp, xc, yc, Xpl, Xpr, Xil, Xir, casechar;
    void* projtmpbuf;
    int *anMaxsR_ = NULL, *anMaxsL_ = NULL, nMaxR, nMaxL;
    int *anMaxsR, *anMaxsL;
    
    // get radii of analyzed ring
    r = 10;
    //R = xs/2-4;
    R = xs / 2;
    sz = 0;
    // evaluate radius
    //    res = IPL_PROJ_FindHoughDonatorProjection5(
    
    // check more args
    if ((psPI == NULL) || (psII == NULL) || (psCI == NULL))
    {
        res = ERROR_NULL_POINTER;
        return res;
    }
    anMaxsR_ = (int*)malloc(((xs + 2) * 2) * sizeof(anMaxsR_[0]));
    anMaxsL_ = anMaxsR_ + xs + 2;
    anMaxsR = anMaxsR_ + 1;
    anMaxsL = anMaxsL_ + 1;
    // normalise center
    xc = psCI->xc;
    yc = psCI->yc;
    if (mode&IVIR_ARRPOXIR_PRPRPR)
    {
        IVIR_ProjPreProcess(projR, R);
        IVIR_ProjPreProcess(projL, R);
    }
    // uint8* medline = (uint8*)malloc(R * sizeof(uint8));
    // printf("PrPu: Calc median\n");
//     res = IPL_HIST_CalcMedianInLine(medline, im, xs, ys, yc, BLUR_HIST_WND);
//     // printf("PrPu: Calc median\n");
//     if (res != ERROR_OK)
//     {
//         free(medline);
//         free(anMaxsR_);
//         return res;
//     }
    
    // list maxima
    nMaxL = ownIVIR_ListMaxs(anMaxsL, projL, r, R);
    nMaxR = ownIVIR_ListMaxs(anMaxsR, projR, r, R);
    // 
    casechar = 0;
    if (mode&IVIR_ARRPOXIR_LEFEYE)
        casechar = 'L';
    if (mode&IVIR_ARRPOXIR_RIGEYE)
        casechar = 'R';
    ownIVIR_SelectPupir(
        (char*)(&casechar), &Xpl, &Xpr, &(psPI->q1), &Xil, &Xir, &(psII->q1),
        projR, anMaxsR, nMaxR, projL, anMaxsL, nMaxL, R, im, xs, ys, xc, yc, mode,
        br_size, br_spread);
    // printf("PrPu: Select PupIr: %d %d %c\n", Xpr, Xpl, casechar);

    /*
    if (imv)
    {
        char nama[FILENAME_MAX];
        DBGL_PUMP_Make3Bpp(imv, xs * 3, im, xs, xs, ys);
        // center, search bound, true pupil, true iris
        DBGL_DRAW_CircleInRGB(imv, xs * 3, xs, ys, xc, yc, 10, 0x00ff00);
        DBGL_DRAW_CircleInRGB(imv, xs * 3, xs, ys, xc, yc, R, 0x00ff00);
        DBGL_DRAW_CircleInRGB(imv, xs * 3, xs, ys, psPI->xc, psPI->yc, psPI->r, 0x003f00);
        DBGL_DRAW_CircleInRGB(imv, xs * 3, xs, ys, psII->xc, psII->yc, psII->r, 0x003f00);
        // detected pupil, detected iris
        DBGL_DRAW_CircleInRGB(imv, xs * 3, xs, ys, xc + (Xpr - Xpl) / 2, yc, (Xpr + Xpl) / 2, 0xff0000);
        DBGL_DRAW_CircleInRGB(imv, xs * 3, xs, ys, xc + (Xir - Xil) / 2, yc, (Xir + Xil) / 2, 0xff7f00);
        // histograms with maxima and principal maximum
        DBGL_DRAW_HistogramInRGBext(imv, R, 100, xs * 3, projL_b, R, 1);
        DBGL_DRAW_HistogramInRGBext(imv + (xs - R) * 3, R, 100, xs * 3, projR_b, R, 1);
        DBGL_DRAW_StrobesInRGB(imv, R, 100, xs * 3, projL_b, R, anMaxsL_, nMaxL * 2 + 1, 0xff0000);
        DBGL_DRAW_StrobesInRGB(imv + (xs - R) * 3, R, 100, xs * 3, projR_b, R, anMaxsR_, nMaxR * 2 + 1, 0xff0000);
        DBGL_DRAW_HistogramInRGBext_uchar(imv + xs*(ys - 50) * 3, xs, 50, xs * 3, medline, xs, 1);
        if (Xpl >= 0)
            anMaxsL[0] = Xpl;
        if (Xil >= 0)
            anMaxsL[(Xpl >= 0)] = Xil;
        if (Xpr >= 0)
            anMaxsR[0] = Xpr;
        if (Xir >= 0)
            anMaxsR[(Xpr >= 0)] = Xir;
        DBGL_DRAW_StrobesInRGB(imv, R, 100, xs * 3, projL_b, R, anMaxsL, (Xpl >= 0) + (Xil >= 0), 0x0000ff);
        DBGL_DRAW_StrobesInRGB(imv + (xs - R) * 3, R, 100, xs * 3, projR_b, R, anMaxsR, (Xpr >= 0) + (Xir >= 0), 0x0000ff);
        // text report
        sprintf(nama, "%c %d", (char)casechar, psII->q2);
        DBGL_TYPE_TextInRGB(imv, xs, ys, nama, 0, ys / 2, 0xff00ff, -1);
        sprintf(nama, "%s_%02x%s.bmp", nambeg, casechar & 0xff, (psII->q2 <= 10) ? "good" : "bad");
        DBGL_FILE_SaveRGB8Image(imv, xs, xs, ys, nama);
        free(imv);
    }
    */
    if ((Xpr >= 0) && (Xpl >= 0))
    {
        psPI->xc = xc + (Xpr - Xpl) / 2;
        psPI->yc = yc;
        psPI->r = (Xpr + Xpl) / 2;
    }
    else
        psPI->xc = psPI->yc = psPI->r = -1;
    // printf("PrPu: Pupil data: %d %d %d\n", psPI->xc, psPI->yc, psPI->r);
    if ((Xir >= 0) && (Xil >= 0))
    {
        psII->xc = xc + (Xir - Xil) / 2;
        psII->yc = yc;
        psII->r = (Xir + Xil) / 2;
    }
    else
        psII->xc = psII->yc = psII->r = -1;
    // printf("PrPu: Iris data: %d %d %d\n", psII->xc, psII->yc, psII->r);
    psII->q2 = casechar;
    free(anMaxsR_);
    return res;
}


/*
static RESULT_CODE ownIVIR_SelectPupirRE(
    char* casechar,   // OUT: character for diagnose
    int* pXpl,
    int* pXpr,
    int *pQp,
    int* pXil,
    int* pXir,
    int* pQi,
    const int* projR,
    const int* anMaxsR,
    int nMaxR,
    const int* projL,
    const int* anMaxsL,
    int nMaxL,
    int projlen,
    const unsigned char* im,
    int xs,
    int ys,
    int xc,
    int yc,
    int mode)
{
    int Xpl, Xpr, Xil, Xir, i, j, k, l, rp, ri, Nr, Nl;
    int Qb, Qd, Qv, Q, max3r, max3l, zs, maxXpl, maxXpr, maxXil, maxXir;
    int valXpl, valXpr, valXil, valXir, PosMaxR, PosMaxL;
    SBazRadHyp* psbrh = NULL;
    int nbrh_used = 0, nbrh_alloced = 0;

    zs = (xs<ys) ? xs : ys;
    // default - not found
    *pXpl = -1;
    *pXpr = -1;
    *pXil = -1;
    *pXir = -1;
    *casechar = '-';
    *pQp = *pQi = -1;
    // init positions
    maxXpl = maxXpr = maxXil = maxXir = -1;
    if ((nMaxR<1) || (nMaxL<1))
    { // no peaks found at least at one side (virtually impossible)
        *casechar = 'N';
        goto quita;
    }
    // find position of max in left and right
    PosMaxR = anMaxsR[0];
    for (i = 1; i<nMaxR; i++)
        if (projR[PosMaxR]<projR[anMaxsR[2 * i]])
            PosMaxR = anMaxsR[2 * i];
    PosMaxL = anMaxsL[0];
    for (i = 1; i<nMaxL; i++)
        if (projL[PosMaxL]<projL[anMaxsL[2 * i]])
            PosMaxL = anMaxsL[2 * i];
    if ((nMaxR<2) || (nMaxL<2))
        // only one peak found (very low probability of this)
        goto bongo;
    // two or more peaks in each side
    for (i = 0; i<nMaxR - 1; i++)
    {
        Xpr = anMaxsR[2 * i];
        for (j = i + 1; j<nMaxR; j++)
        {
            Xir = anMaxsR[2 * j];
            for (k = 0; k<nMaxL - 1; k++)
            {
                Xpl = anMaxsL[2 * k];
                for (l = k + 1; l<nMaxL; l++)
                {
                    Xil = anMaxsL[2 * l];
                    // gupta cheqa
                    if (((projL[PosMaxL]>4 * projL[Xpl]) && (projL[PosMaxL]>4 * projL[Xil])) ||
                        ((projR[PosMaxR]>4 * projR[Xpr]) && (projR[PosMaxR]>4 * projR[Xir])))
                        continue;   // valuable max is out of consideration
                                    // absolute max should be at least one of the borders
                    rp = (Xpl + Xpr) / 2;
                    ri = (Xil + Xir) / 2;
                    // filters
                    if (rp * 7<ri)  // mia120427
                        continue;   // pupil should be at least 1/7 of iris
                    if (rp * 4>ri * 3)
                        continue;   // pupil should be at most 3/4 of iris
                    if (2 * (Xil - Xpl)<(Xir - Xpr))
                        continue;   // iris sides should not differ strongly
                    if (2 * (Xir - Xpr)<(Xil - Xpl))
                        continue;   // iris sides should not differ strongly
                    if ((ri<zs / 10) || (ri>zs / 2) || (rp<5) || (rp>zs / 3))
                        continue;
                    // calculate peak height ratios
                    valXpl = projL[Xpl];
                    valXpr = projR[Xpr];
                    valXil = projL[Xil];
                    valXir = projR[Xir];
                    //== qualities
                    // - pupil/iris decentration
                    if (Xir - Xpr<Xil - Xpl)
                        Qd = 100 * (Xir - Xpr) / (Xil - Xpl);
                    else
                        Qd = 100 * (Xil - Xpl) / (Xir - Xpr);
                    if (Qd<55)
                        Qd = 0;
                    else
                        if (Qd<67)
                            Qd = (Qd - 55)*(100 - 0) / (67 - 55) + 0;
                        else
                            Qd = 100;
                    // - center/iris decentration
                    Qb = 50 * (Xir - Xil) / (Xil + Xir);
                    if (Qb<0)
                        Qb = -Qb;
                    Qb = 100 - Qb;
                    // -value
                    Qv = 50 * (valXpr + valXpl + valXir + valXil) / (projL[PosMaxL] + projR[PosMaxR]);
                    Q = Qv*Qd*Qb / 100 + 1;
                    if (nbrh_used == nbrh_alloced)
                        psbrh = (SBazRadHyp*)realloc(psbrh, (nbrh_alloced += 1024) * sizeof(psbrh[0]));
                    psbrh[nbrh_used].q = Q;
                    psbrh[nbrh_used].xpl = Xpl;
                    psbrh[nbrh_used].xpr = Xpr;
                    psbrh[nbrh_used].xil = Xil;
                    psbrh[nbrh_used].xir = Xir;
                    psbrh[nbrh_used].vpl = valXpl;
                    psbrh[nbrh_used].vpr = valXpr;
                    psbrh[nbrh_used].vil = valXil;
                    psbrh[nbrh_used].vir = valXir;
                    psbrh[nbrh_used].casechar = '0';
                    nbrh_used++;
                }
            }
        }
    }
bongo:;
    // processing for outbounds
    if ((xc<projlen) || (xc + projlen >= xs))
    {
        ownIVIR_GetMaxSubMax_RE(&maxXpl, &maxXil, projL, anMaxsL, nMaxL, &max3l);
        ownIVIR_GetMaxSubMax_RE(&maxXpr, &maxXir, projR, anMaxsR, nMaxR, &max3r);
        // convert index to position
        maxXpl = (maxXpl >= 0) ? anMaxsL[2 * maxXpl] : 0;
        maxXpr = (maxXpr >= 0) ? anMaxsR[2 * maxXpr] : 0;
        maxXil = (maxXil >= 0) ? anMaxsL[2 * maxXil] : 0;
        maxXir = (maxXir >= 0) ? anMaxsR[2 * maxXir] : 0;
        max3l = (max3l >= 0) ? anMaxsL[2 * max3l] : 0;
        max3r = (max3r >= 0) ? anMaxsR[2 * max3r] : 0;
        if (xc<projlen)
        { // left side may be occluded
            if ((max3r<0) ||
                ((2 * projR[maxXpr]>3 * projR[max3r]) && (2 * projR[maxXir]>3 * projR[max3r])))
            { // two right side maxima are good
                valXpr = projR[Xpr = maxXpr];
                valXir = projR[Xir = maxXir];
                valXpl = projL[Xpl = PosMaxL];
                valXil = valXir / 2;
                Xil = Xpl + (Xir - Xpr);
                rp = (Xpl + Xpr) / 2;
                ri = (Xil + Xir) / 2;
                // filters
                if (rp * 7<ri)  //mia120427
                    goto sorry_left;   // pupil should be at least 1/7 of iris
                if (rp * 4>ri * 3)
                    goto sorry_left;   // pupil should be at most 3/4 of iris
                if (2 * (Xil - Xpl)<(Xir - Xpr))
                    goto sorry_left;   // iris sides should not differ strongly
                if (2 * (Xir - Xpr)<(Xil - Xpl))
                    goto sorry_left;   // iris sides should not differ strongly
                if ((ri<zs / 10) || (ri>zs / 2) || (rp<5) || (rp>zs / 3))
                    goto sorry_left;
                if ((Xil <= xc) || (Xil>projlen))
                    goto sorry_left;
                // qualities
                Qb = 50 * (Xir - Xil) / (Xil + Xir);
                if (Qb<0)
                    Qb = -Qb;
                Qb = 100 - Qb;
                Qv = 50 * (valXpr + valXpl + valXir + valXil) / (projL[PosMaxL] + projR[PosMaxR]);
                Q = Qv*Qb * 100 / 100 + 1;
                // make hyp
                if (nbrh_used == nbrh_alloced)
                    psbrh = (SBazRadHyp*)realloc(psbrh, (nbrh_alloced += 64) * sizeof(psbrh[0]));
                psbrh[nbrh_used].q = Q;
                psbrh[nbrh_used].xpl = Xpl;
                psbrh[nbrh_used].xpr = Xpr;
                psbrh[nbrh_used].xil = Xil;
                psbrh[nbrh_used].xir = Xir;
                psbrh[nbrh_used].vpl = valXpl;
                psbrh[nbrh_used].vpr = valXpr;
                psbrh[nbrh_used].vil = valXil;
                psbrh[nbrh_used].vir = valXir;
                psbrh[nbrh_used].casechar = '5';
                nbrh_used++;
            }
        sorry_left:;
        }
        if (xc + projlen >= xs)
        { // left side may be occluded
            if ((max3l<0) ||
                ((2 * projL[maxXpl]>3 * projL[max3l]) && (2 * projL[maxXil]>3 * projL[max3l])))
            { // two right side maxima are good
                valXpl = projL[Xpl = maxXpl];
                valXil = projL[Xil = maxXil];
                valXpr = projR[Xpr = PosMaxR];
                valXir = valXil / 2;
                Xir = Xpr + (Xil - Xpl);
                rp = (Xpl + Xpr) / 2;
                ri = (Xil + Xir) / 2;
                // filters
                if (rp * 7<ri)  //mia120427
                    goto sorry_right;   // pupil should be at least 1/7 of iris
                if (rp * 4>ri * 3)
                    goto sorry_right;   // pupil should be at most 3/4 of iris
                if (2 * (Xil - Xpl)<(Xir - Xpr))
                    goto sorry_right;   // iris sides should not differ strongly
                if (2 * (Xir - Xpr)<(Xil - Xpl))
                    goto sorry_right;   // iris sides should not differ strongly
                if ((ri<zs / 10) || (ri>zs / 2) || (rp<5) || (rp>zs / 3))
                    goto sorry_right;
                // qualities
                if ((Xir + xc<xs) || (Xir>projlen))
                    goto sorry_right;
                Qb = 50 * (Xir - Xil) / (Xil + Xir);
                if (Qb<0)
                    Qb = -Qb;
                Qb = 100 - Qb;
                Qv = 50 * (valXpr + valXpl + valXir + valXil) / (projL[PosMaxL] + projR[PosMaxR]);
                Q = Qv*Qb * 100 / 100 + 1;
                // make hyp
                if (nbrh_used == nbrh_alloced)
                    psbrh = (SBazRadHyp*)realloc(psbrh, (nbrh_alloced += 64) * sizeof(psbrh[0]));
                psbrh[nbrh_used].q = Q;
                psbrh[nbrh_used].xpl = Xpl;
                psbrh[nbrh_used].xpr = Xpr;
                psbrh[nbrh_used].xil = Xil;
                psbrh[nbrh_used].xir = Xir;
                psbrh[nbrh_used].vpl = valXpl;
                psbrh[nbrh_used].vpr = valXpr;
                psbrh[nbrh_used].vil = valXil;
                psbrh[nbrh_used].vir = valXir;
                psbrh[nbrh_used].casechar = '5';
                nbrh_used++;
            }
        }
    sorry_right:;
    }
    // hypotheses are enumerated. Now choose. 
    if (!nbrh_used)
        // no hypothesis at all - nothing found
        goto quita;
    if (nbrh_used == 1)
    { // only one hypothesis
        *casechar = psbrh[0].casechar + 4;
        *pXpl = psbrh[0].xpl;
        *pXpr = psbrh[0].xpr;
        *pXil = psbrh[0].xil;
        *pXir = psbrh[0].xir;
        *pQp = *pQi = 50;
        goto quita;
    }
    qsort(psbrh, nbrh_used, sizeof(psbrh[0]), ownIVIR_sortSBRH_RE);
    // calc ratio of pupil and iris peak heights
    Nr = psbrh[0].vpr * 1024 / (psbrh[0].vir + 1);
    Nl = psbrh[0].vpl * 1024 / (psbrh[0].vil + 1);
    if ((Nr>1024 / 3) && (Nr<1024 * 3) && (Nl>1024 / 3) && (Nl<1024 * 3) &&
        (1000 * psbrh[1].q <= 900 * psbrh[0].q))
    { // almost warrantied
        *casechar = psbrh[0].casechar;//   '*';
        *pXpl = psbrh[0].xpl;
        *pXpr = psbrh[0].xpr;
        *pXil = psbrh[0].xil;
        *pXir = psbrh[0].xir;
        Q = 100;
        *pQp = *pQi = Q;
        goto quita;
    }
    if ((Nr>1024 / 8) && (Nr<1024 * 8) && (Nl>1024 / 8) && (Nl<1024 * 8) &&
        (1000 * psbrh[1].q <= 915 * psbrh[0].q))
    { // bit worse
        *casechar = psbrh[0].casechar + 1;//'+';
        *pXpl = psbrh[0].xpl;
        *pXpr = psbrh[0].xpr;
        *pXil = psbrh[0].xil;
        *pXir = psbrh[0].xir;
        Q = 95;
        *pQp = *pQi = Q;
        goto quita;
    }
    if ((Nr>1024 / 2) && (Nr<1024 * 2) && (Nl>1024 / 2) && (Nl<1024 * 2) &&
        (1000 * psbrh[1].q <= 930 * psbrh[0].q))
    { // worse
        *casechar = psbrh[0].casechar + 2;//'~';
        *pXpl = psbrh[0].xpl;
        *pXpr = psbrh[0].xpr;
        *pXil = psbrh[0].xil;
        *pXir = psbrh[0].xir;
        Q = 90;
        *pQp = *pQi = Q;
        goto quita;
    }
    if ((Nr>1024 / 4) && (Nr<1024 * 4) && (Nl>1024 / 4) && (Nl<1024 * 4) &&
        (1000 * psbrh[1].q <= 930 * psbrh[0].q))
    { // more worse
        *casechar = psbrh[0].casechar + 3;//'!';
        *pXpl = psbrh[0].xpl;
        *pXpr = psbrh[0].xpr;
        *pXil = psbrh[0].xil;
        *pXir = psbrh[0].xir;
        Q = 85;
        *pQp = *pQi = Q;
        goto quita;
    }
    if ((Nr>1024 / 3) && (Nr<1024 * 3) && (Nl>1024 / 3) && (Nl<1024 * 3) &&
        ((psbrh[1].xpl == psbrh[0].xpl) && (psbrh[1].xpr == psbrh[0].xpr) &&
        ((psbrh[1].xil == psbrh[0].xil) || (psbrh[1].xir == psbrh[0].xir))))
    { // max and submax pupils match
        *casechar = 'A';
        *pXpl = psbrh[0].xpl;
        *pXpr = psbrh[0].xpr;
        Q = 75;
        *pQp = Q;
        *pQi = 0;
        goto quita;
    }
    if ((Nr>1024 / 7) && (Nr<1024 * 7) && (Nl>1024 / 7) && (Nl<1024 * 7) &&
        ((psbrh[1].xil == psbrh[0].xil) && (psbrh[1].xir == psbrh[0].xir) &&
        ((psbrh[1].xpl == psbrh[0].xpl) || (psbrh[1].xpr == psbrh[0].xpr))))
    { // max and submax irises match
        *casechar = 'B';
        *pXil = psbrh[0].xil;
        *pXir = psbrh[0].xir;
        Q = 80;
        *pQi = Q;
        *pQp = 0;
        goto quita;
    }
    // just make a silly guess
    ownIVIR_GetMaxSubMax_RE(&maxXpl, &maxXil, projL, anMaxsL, nMaxL, &max3l);
    ownIVIR_GetMaxSubMax_RE(&maxXpr, &maxXir, projR, anMaxsR, nMaxR, &max3r);
    if (maxXpl < 0)
        maxXpl = 0;
    if (maxXpr < 0)
        maxXpr = 0;
    if ((ownIVIR_difil_RE(anMaxsL[2 * maxXil], anMaxsR[2 * maxXir])<150) &&
        (3 * projL[anMaxsL[2 * maxXil]]>2 * projL[PosMaxL]) &&
        (3 * projR[anMaxsR[2 * maxXir]]>2 * projR[PosMaxR]))
    { // second maxima positions collimate, both maxima values big => iris
        if (anMaxsL[2 * maxXil] + anMaxsR[2 * maxXir]>75 * 2)
        {
            *casechar = 'I';
            *pXil = anMaxsL[2 * maxXil];
            *pXir = anMaxsR[2 * maxXir];
            Q = (mode&IVIR_ARRPOXIR_PRPRPR) ? 45 : 83;
            *pQi = Q;
            *pQp = 0;
        }
        else
        {
            *casechar = 'J';
            *pXpl = anMaxsL[2 * maxXil];
            *pXpr = anMaxsR[2 * maxXir];
            //      Q = 100-400*submaxQ/maxQ;
            //    if (Q<=0)
            Q = (mode&IVIR_ARRPOXIR_PRPRPR) ? 10 : 22;
            *pQi = 0;
            *pQp = Q;
        }
        goto quita;
    }
    if ((ownIVIR_difil_RE(anMaxsL[2 * maxXpl], anMaxsR[2 * maxXpr])<150) &&
        (4 * projL[anMaxsL[2 * maxXpl]] >= 2 * projL[PosMaxL]) &&
        (4 * projR[anMaxsR[2 * maxXpr]] >= 2 * projR[PosMaxR]))
    { // first maxima positions collimate, both maxima values big =>
      // this is one circle, but pupil or iris?
        int inL, outL, inR, outR, HL, HR;

        for (inL = i = 0; i<maxXpl; i++)
            inL += projL[anMaxsL[i * 2]];
        for (inR = i = 0; i<maxXpr; i++)
            inR += projR[anMaxsR[i * 2]];
        for (outL = 0, i = maxXpl + 1; i<nMaxL; i++)
            outL += projL[anMaxsL[i * 2]];
        for (outR = 0, i = maxXpr + 1; i<nMaxR; i++)
            outR += projR[anMaxsR[i * 2]];
        HL = projL[anMaxsL[maxXpl * 2]];
        HR = projR[anMaxsR[maxXpr * 2]];
        if ((4 * (outL + outR) <= 5 * (HL + HR)) && (inR + inL >= (HL + HR) / 5))
            //      if (0)
        {
            *casechar = 'X';
            *pXil = anMaxsL[2 * maxXpl];
            *pXir = anMaxsR[2 * maxXpr];
            Q = 45;
            *pQi = Q;
            *pQp = 0;
        }
        else
        { // sum of maxima inside is too small, or outsize too big => pupil
            *pXpl = anMaxsL[2 * maxXpl] - 2;
            *pXpr = anMaxsR[2 * maxXpr] - 2;
            if (4 * (outL + outR) <= 5 * (HL + HR))
            {
                *casechar = 'Z';
                Q = (mode&IVIR_ARRPOXIR_PRPRPR) ? 55 : 70;
            }
            else
            {
                *casechar = 'z';
                Q = (mode&IVIR_ARRPOXIR_PRPRPR) ? 50 : 25;
            }
            *pQp = Q;
            *pQi = 0;
        }
        goto quita;
    }
    if ((ownIVIR_difil_RE(anMaxsL[2 * maxXpl], anMaxsR[2 * maxXir])<100) &&
        (3 * projL[anMaxsL[2 * maxXpl]]>2 * projL[PosMaxL]) &&
        (3 * projR[anMaxsR[2 * maxXir]]>2 * projR[PosMaxR]))
    { // consider as iris
        j = anMaxsL[2 * maxXpl];
        i = anMaxsR[2 * maxXir];
        // !!! size hack - decision depends on image size !!!
        if (((i<projlen / 4) && (j<projlen / 4)) || (i + j<75))
        {
            Q = (mode&IVIR_ARRPOXIR_PRPRPR) ? 35 : 20;
            *casechar = 'U';
            *pXpr = i - 2;
            *pXpl = j - 2;
            *pQp = Q;
            *pQi = 0;
        }
        else
        {
            Q = (mode&IVIR_ARRPOXIR_PRPRPR) ? 30 : 20;
            *casechar = 'u';
            *pXir = i;
            *pXil = j;
            *pQi = Q;
            *pQp = 0;
        }
        goto quita;
    }
    if ((ownIVIR_difil_RE(anMaxsL[2 * maxXil], anMaxsR[2 * maxXpr])<100) &&
        (3 * projL[anMaxsL[2 * maxXil]]>2 * projL[PosMaxL]) &&
        (3 * projR[anMaxsR[2 * maxXpr]]>2 * projR[PosMaxR]))
    { // consider as iris
        j = anMaxsL[2 * maxXil];
        i = anMaxsR[2 * maxXpr];
        // !!! size hack - decision depends on image size !!!
        if ((i<projlen / 4) && (j<projlen / 4))
        {
            Q = (mode&IVIR_ARRPOXIR_PRPRPR) ? 40 : 10;
            *casechar = 'V';
            *pXpr = i - 2;
            *pXpl = j - 2;
            *pQp = Q;
            *pQi = 0;
        }
        else
        {
            Q = (mode&IVIR_ARRPOXIR_PRPRPR) ? 25 : 10;
            *casechar = 'v';
            *pXir = i;
            *pXil = j;
            *pQi = Q;
            *pQp = 0;
        }
        goto quita;
    }
    i = anMaxsR[2 * ((projR[anMaxsR[2 * maxXpr]]>projR[anMaxsR[2 * maxXir]]) ? maxXpr : maxXir)];
    j = anMaxsL[2 * ((projL[anMaxsL[2 * maxXpl]]>projL[anMaxsL[2 * maxXil]]) ? maxXpl : maxXil)];
    // !!! size hack - decision depends on image size !!!
    if ((i<2 * 3 * projlen / (5 * 4)) && (j<2 * 3 * projlen / (5 * 4)))
    {
        *casechar = 'g';
        *pXpr = i - 2;
        *pXpl = j - 2;
        *pQp = (mode&IVIR_ARRPOXIR_PRPRPR) ? 20 : 15;
        *pQi = 0;
    }
    else
    {
        *casechar = 'G';
        *pXir = i;
        *pXil = j;
        *pQi = (mode&IVIR_ARRPOXIR_PRPRPR) ? 15 : 2;
        *pQp = 0;
    }
quita:
    if (psbrh)
        free(psbrh);
    return ERROR_OK;
}*/

/*
// calculate left side, right side and total circular projection 
// in a concentric ring
MIA_RESULT_CODE IPL_PROJ_FindHoughDonatorProjection3(
  int* pnProjT,             // OUT: circular projection - total
  int* pnProjR,             // OUT: circular projection - right side
  int* pnProjL,             // OUT: circular projection - left side
  const unsigned char* im,  // IN:  image
  int xs,                   // IN:  image size
  int ys,                   // 
  int xc,                   // IN:  ring center
  int yc,                   //
  int r,                    // IN:  inner radius
  int* pR,                  // IN/OUT: outer radius proposed/actually used
  int nMaxRadOfDence,       // IN:  maximum radius for dence processing. 0 - dence always
  void* buf,                // IN:  external buffer
  int* buflen)              // IN/OUT: allocated/used bytes
{
  int i,j,gx,gy,y,G,L,S,*numpoR,*numpoL,*numpoT,currad,R;
  int j_inn,j_out,Y,D,RR_ii,rr_ii;
  const unsigned char* pim;

    // check arguments
    if ((buflen==NULL)||(pR==NULL))
      return ERR_GEN_NULLPOINTER;
    R = *pR;
    if (R<10)
      return ERR_GEN_SIZE_NOT_MATCH;
    // calc size
    i = (R+1)*3*sizeof(numpoT[0]);//+sizeof(short)*xs*ys+sizeof(int)*(xs+ys);
    if (buf==NULL)
    {
      *buflen = i;
      return ERR_OK;
    }
    if (*buflen<i)
    {
      *buflen = i;
      return ERR_GEN_NOMEMORY;
    }
    // more check arguments
    if ((pnProjT==NULL)||(pnProjR==NULL)||(pnProjL==NULL)||(im==NULL))
      return ERR_GEN_NULLPOINTER;
    if ((((ptr_t)im)&3)||(((ptr_t)buf)&3)||(xs&3)||
        (((ptr_t)pnProjT)&3)||(((ptr_t)pnProjR)&3)||(((ptr_t)pnProjL)&3))
      return ERR_GEN_BAD_ALIGNMENT;
    if ((r+5>R)||(r<5))
      return ERR_GEN_SIZE_NOT_MATCH;
    if ((xs<2*r)||(ys<2*r))
      return ERR_GEN_NO_DATA;
    if ((xc-r/2<0)||(xc+r/2>=xs)||(yc-r/2<0)||(yc+r/2>=ys))
      return ERR_GEN_INVALID_PARAMETER;
    // verify outer radius - bound by maximum span from center to border
    i = xs-1-xc;
    if (i<xc)
      i = xc;
    if (R>i)
      R = i;
    // verify outer radius - bound by image size
    if (2*R>ys)
      R = ys/2;
    if (2*R>xs)
      R = xs/2;
    // return corrected outer radius
    *pR = R;
    // allocate
    numpoT = (int*)buf;
    numpoR = numpoT+(R+1);
    numpoL = numpoR+(R+1);
    // clear junk
    MIA_memset(pnProjT,0,sizeof(pnProjT[0])*(1+R));
    MIA_memset(pnProjR,0,sizeof(pnProjR[0])*(1+R));
    MIA_memset(pnProjL,0,sizeof(pnProjL[0])*(1+R));
    MIA_memset(numpoT,0,3*sizeof(numpoT[0])*(1+R));
    // along y
    if (!nMaxRadOfDence)
      D = 0x10000;
    else
    {
      D = (r*0x10000)/nMaxRadOfDence;
      if (D<0x10000)
        D = 0x10000;
    }
    for (Y=-R*0x10000;Y<=R*0x10000;Y+=D)
    {
      i = (Y+0x8000)/0x10000;
      y = yc+i;
      if ((y<1)||(y>=ys-1))
        continue;
      // left side
      j_out = RR_ii = MIA_INT_sqrt(R*R-i*i);
      if (xc-j_out<1)
        j_out = xc-1;
      j_inn = (i>0)?i:(-i);
      if (j_inn<r)
      {
        rr_ii = MIA_INT_sqrt(r*r-i*i);
        if (j_inn<rr_ii)
          j_inn = rr_ii;
      }
      L = i*i+j_out*j_out;
      currad = MIA_INT_sqrt(L);
      pim = im+y*xs+xc-j_out;
      for (j=-j_out;j<=-j_inn;j++)
      {
        // count number of points in circle
        numpoT[currad]++;
        numpoL[currad]++;
        // straight gradient
        gx =+pim[xs+1]+pim[ 1]+pim[-xs+1]
            -pim[xs-1]-pim[-1]-pim[-xs-1];
        if (gx<-4)
        {
          gy =+pim[xs+1]+pim[xs]+pim[xs-1]
              -pim[-xs+1]-pim[-xs]-pim[-xs-1];
          G = gx*gx+gy*gy;
          if (G>36*2)
          {
            // inner product
            S = gx*j+gy*i;
            //
            if ((S>0)&&(((int40)(S*2))*S>((int40)L)*G))
            {
              pnProjT[currad]++;
              pnProjL[currad]++;
            }
          }
        }
        pim++;
        L -= (-j*2-1);
        if (L<currad*currad)
          currad--;
      }
      // right side
      j_out = RR_ii;
      if (xc+j_out>xs-2)
        j_out = xs-2-xc;
      L = i*i+j_out*j_out;
      currad = MIA_INT_sqrt(L);
      pim = im+y*xs+xc+j_out;
      for (j=j_out;j>=j_inn;j--)
      {
        // count number of points in circle
        numpoT[currad]++;
        numpoR[currad]++;
        // straight gradient
        gx =+pim[xs+1]+pim[ 1]+pim[-xs+1]
            -pim[xs-1]-pim[-1]-pim[-xs-1];
        if (gx>4)
        {
          gy =+pim[xs+1]+pim[xs]+pim[xs-1]
              -pim[-xs+1]-pim[-xs]-pim[-xs-1];
          G = gx*gx+gy*gy;
          if (G>36*2)
          {
            // inner product
            S = gx*j+gy*i;
            //
            if ((S>0)&&(((int40)(S*2))*S>((int40)L)*G))
            {
              pnProjT[currad]++;
              pnProjR[currad]++;
            }
          }
        }
        pim--;
        L -= (j*2-1);
        if (L<currad*currad)
          currad--;
      }
    }
    // normalize
    for (i=0;i<=R;i++)
    {
      pnProjT[i] = (pnProjT[i]*1024)/(numpoT[i]+1);
      pnProjR[i] = (pnProjR[i]*1024)/(numpoR[i]+1);
      pnProjL[i] = (pnProjL[i]*1024)/(numpoL[i]+1);
    }
    return ERR_OK;
}

#define BLWIN_I2P 5

// calculate four quadrant projections (for pupil-from-iris)
MIA_RESULT_CODE IPL_PROJ_FindHoughDonatorProjection4_v2(
  int* pnProjR,             // OUT: circular projection - right side
  int* pnProjL,             // OUT: circular projection - left side
  int* pnProjT,             // OUT: circular projection - top side
  int* pnProjB,             // OUT: circular projection - bottom side
  const unsigned char* im,  // IN:  image
  int xs,                   // IN:  image size
  int ys,                   // 
  int xc,                   // IN:  center
  int yc,                   //
  int R,                    // IN:  outer radius
  void* buf,                // IN:  temporary buffer
  int* plen,                // IN/OUT: temporary buffer allocated/used
  const char* nambeg)
{
  MIA_RESULT_CODE res=ERR_OK;
  int *pnProjR_n,*pnProjL_n,*pnProjT_n,*pnProjB_n,*pnProjX,*pnProjX_n;
  int x,y,w,h,cx,cy,gx,gy,cIdx,sz;
  int64 G,L,P,ll;
  uint8 *dbgim=NULL;
//  char nama[FILENAME_MAX];
  double S,Sg,Sgg;
  int T;

    // check arguments
    if (plen==NULL)
      return ERR_GEN_NULLPOINTER;
    if (R<=0)
      return ERR_GEN_NO_DATA;
    // calculate size
    sz = R*6*sizeof(pnProjR_n[0]);
    if (buf==NULL)
    {
      *plen = sz;
      return ERR_OK;
    }
    if (*plen<sz)
      return ERR_GEN_INSUFFICIENT_BUFFER;
    // check more arguments
    if ((pnProjR==NULL)||(pnProjL==NULL)||(pnProjT==NULL)||(pnProjB==NULL)||(im==NULL))
      return ERR_GEN_NULLPOINTER;
    if (((ptr_t)buf)&3)
      return ERR_GEN_BAD_ALIGNMENT;
    if ((xs<2*R)||(ys<2*R))
      return ERR_GEN_BAD_DATA;
    if ((xc-R/10<0)||(xc+R/10>=xs)||(yc-R/10<0)||(yc+R/10>=ys))
      return ERR_GEN_INVALID_PARAMETER; // too close/out of border
    // allocate
    *plen = sz;
    pnProjR_n = (int*)buf;
    pnProjL_n = pnProjR_n+R;
    pnProjT_n = pnProjL_n+R;
    pnProjB_n = pnProjT_n+R;
    pnProjX   = pnProjB_n+R;
    pnProjX_n = pnProjX+R;
    // clear junk
    MIA_memset(pnProjR,0,sizeof(pnProjR[0])*R);
    MIA_memset(pnProjL,0,sizeof(pnProjL[0])*R);
    MIA_memset(pnProjT,0,sizeof(pnProjR[0])*R);
    MIA_memset(pnProjB,0,sizeof(pnProjL[0])*R);
    MIA_memset(pnProjR_n,0,sizeof(pnProjR_n[0])*R*4);
    // obtain dimensions, remember to leave 1-pixel gap for Sobel mask
    x = xc-R;
    y = yc-R;
    w = 2*R+1;
    h = 2*R+1;
    if (x<1)
    {
      w += x-1;
      x = 1;
    }
    if (y<1)
    {
      h += y-1;
      y = 1;
    }
    if (x+w>=xs-1)
      w = xs-x-2;
    if (y+h>=ys-1)
      h = ys-y-2;
    // allocate debug image
    if (nambeg)
    {
      dbgim = (uint8*)malloc(w*h);
      memset(dbgim,0,w*h);
    }
    // calculate gradient threshold
    S = Sg = Sgg = 0.;
    for (cy=y;cy<y+h;cy++)
    {
      for (cx=x;cx<x+w;cx++)
      {
        L = (cx-xc)*(cx-xc)+(cy-yc)*(cy-yc);
        ll = (int)(sqrt((double)L)+.5);
        if ((ll<R/10)||(ll>3*R/4))
          continue;
        gx = im[(cy-1)*xs+(cx+1)]+2*im[(cy  )*xs+(cx+1)]+im[(cy+1)*xs+(cx+1)]
            -im[(cy-1)*xs+(cx-1)]-2*im[(cy  )*xs+(cx-1)]-im[(cy+1)*xs+(cx-1)];
        gy = im[(cy+1)*xs+(cx+1)]+2*im[(cy+1)*xs+(cx  )]+im[(cy+1)*xs+(cx-1)]
            -im[(cy-1)*xs+(cx+1)]-2*im[(cy-1)*xs+(cx  )]-im[(cy-1)*xs+(cx-1)];
        G = gx*gx+gy*gy;
        if (G<100)
          continue;
        P = (cx-xc)*gx+(cy-yc)*gy;
        if (P*P*10<L*G*9)
          continue;
        S++;
        Sg += sqrt((double)G);
        Sgg += G;
      }
    }
    if (S<=1)
    {
      res = ERR_GEN_NO_DATA;
      goto quitti;
    }
    Sg  /= S;
    Sgg /= S;
    //T = (int)(Sg+sqrt(Sgg-Sg*Sg)+.5);
    //T = (int)(Sg+.5);
    T = (int)(Sg-sqrt(Sgg-Sg*Sg)+.5);
    // calculate gradient
    for (cy=y;cy<y+h;cy++)
    {
      for (cx=x;cx<x+w;cx++)
      {
        gx = im[(cy-1)*xs+(cx+1)]+2*im[(cy  )*xs+(cx+1)]+im[(cy+1)*xs+(cx+1)]
            -im[(cy-1)*xs+(cx-1)]-2*im[(cy  )*xs+(cx-1)]-im[(cy+1)*xs+(cx-1)];
        gy = im[(cy+1)*xs+(cx+1)]+2*im[(cy+1)*xs+(cx  )]+im[(cy+1)*xs+(cx-1)]
            -im[(cy-1)*xs+(cx+1)]-2*im[(cy-1)*xs+(cx  )]-im[(cy-1)*xs+(cx-1)];
        L = (cx-xc)*(cx-xc)+(cy-yc)*(cy-yc);
        ll = (int)(sqrt((double)L)+.5);
        if ((cx-xc)+(cy-yc)>0)
        { // right or top
          if ((cx-xc)-(cy-yc)>0)
            pnProjR_n[ll]++;
          else
            pnProjT_n[ll]++;
        }
        else
        { // left or bottom
          if ((cx-xc)-(cy-yc)>0)
            pnProjB_n[ll]++;
          else
            pnProjL_n[ll]++;
        }
        G = gx*gx+gy*gy;
        if (G<100)
          continue;
        P = (cx-xc)*gx+(cy-yc)*gy;
        if (P*P*10<L*G*9)
          continue;
        if (sqrt((double)G)<T)
          continue;
        if ((cx-xc)+(cy-yc)>0)
        { // right or top
          if ((cx-xc)-(cy-yc)>0)
            pnProjR[ll]++;
          else
            pnProjT[ll]++;
        }
        else
        { // left or bottom
          if ((cx-xc)-(cy-yc)>0)
            pnProjB[ll]++;
          else
            pnProjL[ll]++;
        }
        if (dbgim)
          dbgim[(cy-y)*w+(cx-x)] = 255;
      }
    }
    // blur & normalize histograms
    IPL_HIST_Blur(pnProjX,pnProjR,R,BLWIN_I2P);
    IPL_HIST_Blur(pnProjX_n,pnProjR_n,R,BLWIN_I2P);
    for (cIdx=0;cIdx<R;cIdx++)
      pnProjR[cIdx] = (pnProjX_n[cIdx])?(1000*pnProjX[cIdx]/pnProjX_n[cIdx]):0;
    for (cIdx=R/10;(pnProjR[cIdx+1]<=pnProjR[cIdx])&&(cIdx<R);cIdx++);
    for (;cIdx>=0;pnProjR[cIdx--]=1);
    for (cIdx=R*3/4;(pnProjR[cIdx-1]<=pnProjR[cIdx])&&(cIdx);cIdx--);
    for (;cIdx<R;pnProjR[cIdx++]=1);
    IPL_HIST_Blur(pnProjX,pnProjL,R,BLWIN_I2P);
    IPL_HIST_Blur(pnProjX_n,pnProjL_n,R,BLWIN_I2P);
    for (cIdx=0;cIdx<R;cIdx++)
      pnProjL[cIdx] = (pnProjX_n[cIdx])?(1000*pnProjX[cIdx]/pnProjX_n[cIdx]):0;
    for (cIdx=R/10;(pnProjL[cIdx+1]<=pnProjL[cIdx])&&(cIdx<R);cIdx++);
    for (;cIdx>=0;pnProjL[cIdx--]=1);
    for (cIdx=R*3/4;(pnProjL[cIdx-1]<=pnProjL[cIdx])&&(cIdx);cIdx--);
    for (;cIdx<R;pnProjL[cIdx++]=1);
    IPL_HIST_Blur(pnProjX,pnProjB,R,BLWIN_I2P);
    IPL_HIST_Blur(pnProjX_n,pnProjB_n,R,BLWIN_I2P);
    for (cIdx=0;cIdx<R;cIdx++)
      pnProjB[cIdx] = (pnProjX_n[cIdx])?(1000*pnProjX[cIdx]/pnProjX_n[cIdx]):0;
    for (cIdx=R/10;(pnProjB[cIdx+1]<=pnProjB[cIdx])&&(cIdx<R);cIdx++);
    for (;cIdx>=0;pnProjB[cIdx--]=1);
    for (cIdx=R*3/4;(pnProjB[cIdx-1]<=pnProjB[cIdx])&&(cIdx);cIdx--);
    for (;cIdx<R;pnProjB[cIdx++]=1);
    IPL_HIST_Blur(pnProjX,pnProjT,R,BLWIN_I2P);
    IPL_HIST_Blur(pnProjX_n,pnProjT_n,R,BLWIN_I2P);
    for (cIdx=0;cIdx<R;cIdx++)
      pnProjT[cIdx] = (pnProjX_n[cIdx])?(1000*pnProjX[cIdx]/pnProjX_n[cIdx]):0;
    for (cIdx=R/10;(pnProjT[cIdx+1]<=pnProjT[cIdx])&&(cIdx<R);cIdx++);
    for (;cIdx>=0;pnProjT[cIdx--]=1);
    for (cIdx=R*3/4;(pnProjT[cIdx-1]<=pnProjT[cIdx])&&(cIdx);cIdx--);
    for (;cIdx<R;pnProjT[cIdx++]=1);
    // free & exit
quitti:;
    if (dbgim)
      free(dbgim);
    return res;
}

// calculate four quadrant projections (for pupil-from-iris)
MIA_RESULT_CODE IPL_PROJ_FindHoughDonatorProjection4(
  int* pnProjR,             // OUT: circular projection - right side
  int* pnProjL,             // OUT: circular projection - left side
  int* pnProjT,             // OUT: circular projection - top side
  int* pnProjB,             // OUT: circular projection - bottom side
  const unsigned char* im,  // IN:  image
  int xs,                   // IN:  image size
  int ys,                   // 
  int xc,                   // IN:  center
  int yc,                   //
  int R,                    // IN:  outer radius
  void* buf,                // IN:  temporary buffer
  int* plen,                // IN/OUT: temporary buffer allocated/used
  const char* nambeg)
{
  MIA_RESULT_CODE res=ERR_OK;
  int *pnProjR_n,*pnProjL_n,*pnProjT_n,*pnProjB_n,*pnProjX,*pnProjX_n;
  int x,y,w,h,cx,cy,gx,gy,cIdx,sz;
  int64 G,L,P,ll;
  uint8 *dbgim=NULL;
//  char nama[FILENAME_MAX];

    // check arguments
    if (plen==NULL)
      return ERR_GEN_NULLPOINTER;
    if (R<=0)
      return ERR_GEN_NO_DATA;
    // calculate size
    sz = R*6*sizeof(pnProjR_n[0]);
    if (buf==NULL)
    {
      *plen = sz;
      return ERR_OK;
    }
    if (*plen<sz)
      return ERR_GEN_INSUFFICIENT_BUFFER;
    // check more arguments
    if ((pnProjR==NULL)||(pnProjL==NULL)||(pnProjT==NULL)||(pnProjB==NULL)||(im==NULL))
      return ERR_GEN_NULLPOINTER;
    if (((ptr_t)buf)&3)
      return ERR_GEN_BAD_ALIGNMENT;
    if ((xs<2*R)||(ys<2*R))
      return ERR_GEN_BAD_DATA;
    if ((xc-R/10<0)||(xc+R/10>=xs)||(yc-R/10<0)||(yc+R/10>=ys))
      return ERR_GEN_INVALID_PARAMETER; // too close/out of border
    // allocate
    *plen = sz;
    pnProjR_n = (int*)buf;
    pnProjL_n = pnProjR_n+R;
    pnProjT_n = pnProjL_n+R;
    pnProjB_n = pnProjT_n+R;
    pnProjX   = pnProjB_n+R;
    pnProjX_n = pnProjX+R;
    // clear junk
    MIA_memset(pnProjR,0,sizeof(pnProjR[0])*R);
    MIA_memset(pnProjL,0,sizeof(pnProjL[0])*R);
    MIA_memset(pnProjT,0,sizeof(pnProjR[0])*R);
    MIA_memset(pnProjB,0,sizeof(pnProjL[0])*R);
    MIA_memset(pnProjR_n,0,sizeof(pnProjR_n[0])*R*4);
    // obtain dimensions, remember to leave 1-pixel gap for Sobel mask
    x = xc-R;
    y = yc-R;
    w = 2*R+1;
    h = 2*R+1;
    if (x<1)
    {
      w += x-1;
      x = 1;
    }
    if (y<1)
    {
      h += y-1;
      y = 1;
    }
    if (x+w>=xs-1)
      w = xs-x-2;
    if (y+h>=ys-1)
      h = ys-y-2;
    // allocate debug image
    if (nambeg)
    {
      dbgim = (uint8*)malloc(w*h);
      memset(dbgim,0,w*h);
    }
    // calculate gradient
    for (cy=y;cy<y+h;cy++)
    {
      for (cx=x;cx<x+w;cx++)
      {
        gx = im[(cy-1)*xs+(cx+1)]+2*im[(cy  )*xs+(cx+1)]+im[(cy+1)*xs+(cx+1)]
            -im[(cy-1)*xs+(cx-1)]-2*im[(cy  )*xs+(cx-1)]-im[(cy+1)*xs+(cx-1)];
        gy = im[(cy+1)*xs+(cx+1)]+2*im[(cy+1)*xs+(cx  )]+im[(cy+1)*xs+(cx-1)]
            -im[(cy-1)*xs+(cx+1)]-2*im[(cy-1)*xs+(cx  )]-im[(cy-1)*xs+(cx-1)];
        L = (cx-xc)*(cx-xc)+(cy-yc)*(cy-yc);
        ll = (int)(sqrt((double)L)+.5);
        if ((cx-xc)+(cy-yc)>0)
        { // right or top
          if ((cx-xc)-(cy-yc)>0)
            pnProjR_n[ll]++;
          else
            pnProjT_n[ll]++;
        }
        else
        { // left or bottom
          if ((cx-xc)-(cy-yc)>0)
            pnProjB_n[ll]++;
          else
            pnProjL_n[ll]++;
        }
        G = gx*gx+gy*gy;
        if (G<200)
          continue;
        P = (cx-xc)*gx+(cy-yc)*gy;
        if (P*P*10<L*G*9)        //if (P*P*12<L*G*11)
          continue;
        if ((cx-xc)+(cy-yc)>0)
        { // right or top
          if ((cx-xc)-(cy-yc)>0)
            pnProjR[ll]++;
          else
            pnProjT[ll]++;
        }
        else
        { // left or bottom
          if ((cx-xc)-(cy-yc)>0)
            pnProjB[ll]++;
          else
            pnProjL[ll]++;
        }
        if (dbgim)
          dbgim[(cy-y)*w+(cx-x)] = 255;
      }
    }
    // blur & normalize histograms
    IPL_HIST_Blur(pnProjX,pnProjR,R,BLWIN_I2P);
    IPL_HIST_Blur(pnProjX_n,pnProjR_n,R,BLWIN_I2P);
    for (cIdx=0;cIdx<R;cIdx++)
      pnProjR[cIdx] = (pnProjX_n[cIdx])?(1000*pnProjX[cIdx]/pnProjX_n[cIdx]):0;
    for (cIdx=R/10;(pnProjR[cIdx+1]<=pnProjR[cIdx])&&(cIdx<R);cIdx++);
    for (;cIdx>=0;pnProjR[cIdx--]=1);
    for (cIdx=R*3/4;(pnProjR[cIdx-1]<=pnProjR[cIdx])&&(cIdx);cIdx--);
    for (;cIdx<R;pnProjR[cIdx++]=1);
    IPL_HIST_Blur(pnProjX,pnProjL,R,BLWIN_I2P);
    IPL_HIST_Blur(pnProjX_n,pnProjL_n,R,BLWIN_I2P);
    for (cIdx=0;cIdx<R;cIdx++)
      pnProjL[cIdx] = (pnProjX_n[cIdx])?(1000*pnProjX[cIdx]/pnProjX_n[cIdx]):0;
    for (cIdx=R/10;(pnProjL[cIdx+1]<=pnProjL[cIdx])&&(cIdx<R);cIdx++);
    for (;cIdx>=0;pnProjL[cIdx--]=1);
    for (cIdx=R*3/4;(pnProjL[cIdx-1]<=pnProjL[cIdx])&&(cIdx);cIdx--);
    for (;cIdx<R;pnProjL[cIdx++]=1);
    IPL_HIST_Blur(pnProjX,pnProjB,R,BLWIN_I2P);
    IPL_HIST_Blur(pnProjX_n,pnProjB_n,R,BLWIN_I2P);
    for (cIdx=0;cIdx<R;cIdx++)
      pnProjB[cIdx] = (pnProjX_n[cIdx])?(1000*pnProjX[cIdx]/pnProjX_n[cIdx]):0;
    for (cIdx=R/10;(pnProjB[cIdx+1]<=pnProjB[cIdx])&&(cIdx<R);cIdx++);
    for (;cIdx>=0;pnProjB[cIdx--]=1);
    for (cIdx=R*3/4;(pnProjB[cIdx-1]<=pnProjB[cIdx])&&(cIdx);cIdx--);
    for (;cIdx<R;pnProjB[cIdx++]=1);
    IPL_HIST_Blur(pnProjX,pnProjT,R,BLWIN_I2P);
    IPL_HIST_Blur(pnProjX_n,pnProjT_n,R,BLWIN_I2P);
    for (cIdx=0;cIdx<R;cIdx++)
      pnProjT[cIdx] = (pnProjX_n[cIdx])?(1000*pnProjX[cIdx]/pnProjX_n[cIdx]):0;
    for (cIdx=R/10;(pnProjT[cIdx+1]<=pnProjT[cIdx])&&(cIdx<R);cIdx++);
    for (;cIdx>=0;pnProjT[cIdx--]=1);
    for (cIdx=R*3/4;(pnProjT[cIdx-1]<=pnProjT[cIdx])&&(cIdx);cIdx--);
    for (;cIdx<R;pnProjT[cIdx++]=1);
    // free & exit
    if (dbgim)
      free(dbgim);
    return res;
}

// calculate four quadrant projections (for pupil-from-iris)
MIA_RESULT_CODE IPL_PROJ_FindHoughDonatorProjection4_v1(
  int* pnProjR,             // OUT: circular projection - right side
  int* pnProjL,             // OUT: circular projection - left side
  int* pnProjT,             // OUT: circular projection - top side
  int* pnProjB,             // OUT: circular projection - bottom side
  const unsigned char* im,  // IN:  image
  int xs,                   // IN:  image size
  int ys,                   // 
  int xc,                   // IN:  center
  int yc,                   //
  int R,                    // IN:  outer radius
  void* buf,                // IN:  temporary buffer
  int* plen,                // IN/OUT: temporary buffer allocated/used
  const char* nambeg)
{
  MIA_RESULT_CODE res=ERR_OK;
  int *pnProjR_s,*pnProjL_s,*pnProjT_s,*pnProjB_s,*pnProjX_s,*pnProjX;
  int x,y,w,h,cx,cy,gx,gy,cIdx,sz;
  int64 G,L,P,ll;
  uint8 *dbgim=NULL;
//  char nama[FILENAME_MAX];

    // check arguments
    if (plen==NULL)
      return ERR_GEN_NULLPOINTER;
    if (R<=0)
      return ERR_GEN_NO_DATA;
    // calculate size
    sz = R*6*sizeof(pnProjR[0]);
    if (buf==NULL)
    {
      *plen = sz;
      return ERR_OK;
    }
    if (*plen<sz)
      return ERR_GEN_INSUFFICIENT_BUFFER;
    // check more arguments
    if ((pnProjR==NULL)||(pnProjL==NULL)||(pnProjT==NULL)||(pnProjB==NULL)||(im==NULL))
      return ERR_GEN_NULLPOINTER;
    if (((ptr_t)buf)&3)
      return ERR_GEN_BAD_ALIGNMENT;
    if ((xs<2*R)||(ys<2*R))
      return ERR_GEN_BAD_DATA;
    if ((xc-R/10<0)||(xc+R/10>=xs)||(yc-R/10<0)||(yc+R/10>=ys))
      return ERR_GEN_INVALID_PARAMETER; // too close/out of border
    // allocate
    *plen = sz;
    pnProjR_s = (int*)buf;
    pnProjL_s = pnProjR_s+R;
    pnProjT_s = pnProjL_s+R;
    pnProjB_s = pnProjT_s+R;
    pnProjX_s = pnProjB_s+R;
    pnProjX = pnProjX_s+R;
    // clear junk
    MIA_memset(pnProjR,0,sizeof(pnProjR_s[0])*R);
    MIA_memset(pnProjL,0,sizeof(pnProjL_s[0])*R);
    MIA_memset(pnProjT,0,sizeof(pnProjR_s[0])*R);
    MIA_memset(pnProjB,0,sizeof(pnProjL_s[0])*R);
    MIA_memset(pnProjR_s,0,sizeof(pnProjR[0])*R*4);
    // obtain dimensions, remember to leave 1-pixel gap for Sobel mask
    x = xc-R;
    y = yc-R;
    w = 2*R+1;
    h = 2*R+1;
    if (x<1)
    {
      w += x-1;
      x = 1;
    }
    if (y<1)
    {
      h += y-1;
      y = 1;
    }
    if (x+w>=xs-1)
      w = xs-x-2;
    if (y+h>=ys-1)
      h = ys-y-2;
    // allocate debug image
    if (nambeg)
    {
      dbgim = (uint8*)malloc(w*h);
      memset(dbgim,0,w*h);
    }
    // calculate gradient
    for (cy=y;cy<y+h;cy++)
    {
      for (cx=x;cx<x+w;cx++)
      {
        gx = im[(cy-1)*xs+(cx+1)]+2*im[(cy  )*xs+(cx+1)]+im[(cy+1)*xs+(cx+1)]
            -im[(cy-1)*xs+(cx-1)]-2*im[(cy  )*xs+(cx-1)]-im[(cy+1)*xs+(cx-1)];
        gy = im[(cy+1)*xs+(cx+1)]+2*im[(cy+1)*xs+(cx  )]+im[(cy+1)*xs+(cx-1)]
            -im[(cy-1)*xs+(cx+1)]-2*im[(cy-1)*xs+(cx  )]-im[(cy-1)*xs+(cx-1)];
        L = (cx-xc)*(cx-xc)+(cy-yc)*(cy-yc);
        ll = (int)(sqrt((double)L)+.5);
        if ((cx-xc)+(cy-yc)>0)
        { // right or top
          if ((cx-xc)-(cy-yc)>0)
            pnProjR[ll]++;
          else
            pnProjT[ll]++;
        }
        else
        { // left or bottom
          if ((cx-xc)-(cy-yc)>0)
            pnProjB[ll]++;
          else
            pnProjL[ll]++;
        }
        G = gx*gx+gy*gy;
        if (G<200)
          continue;
        P = (cx-xc)*gx+(cy-yc)*gy;
        if (P*P*10<L*G*9)        //if (P*P*12<L*G*11)
          continue;
        if ((cx-xc)+(cy-yc)>0)
        { // right or top
          if ((cx-xc)-(cy-yc)>0)
            pnProjR_s[ll]++;
          else
            pnProjT_s[ll]++;
        }
        else
        { // left or bottom
          if ((cx-xc)-(cy-yc)>0)
            pnProjB_s[ll]++;
          else
            pnProjL_s[ll]++;
        }
        if (dbgim)
          dbgim[(cy-y)*w+(cx-x)] = 255;
      }
    }

    // blur & normalize histograms
    IPL_HIST_Blur(pnProjX_s,pnProjR_s,R,BLWIN_I2P);
    IPL_HIST_Blur(pnProjX,pnProjR,R,BLWIN_I2P);
    for (cIdx=0;cIdx<R;cIdx++)
      pnProjR_s[cIdx] = (pnProjX[cIdx])?(1000*pnProjX_s[cIdx]/pnProjX[cIdx]):0;
    for (cIdx=R/10;(pnProjR_s[cIdx+1]<=pnProjR_s[cIdx])&&(cIdx<R);cIdx++);
    for (;cIdx>=0;pnProjR_s[cIdx--]=1);
    for (cIdx=R*3/4;(pnProjR_s[cIdx-1]<=pnProjR_s[cIdx])&&(cIdx);cIdx--);
    for (;cIdx<R;pnProjR_s[cIdx++]=1);
memcpy(pnProjR,pnProjR_s,R*4);
    IPL_HIST_Blur(pnProjX_s,pnProjL_s,R,BLWIN_I2P);
    IPL_HIST_Blur(pnProjX,pnProjL,R,BLWIN_I2P);
    for (cIdx=0;cIdx<R;cIdx++)
      pnProjL_s[cIdx] = (pnProjX[cIdx])?(1000*pnProjX_s[cIdx]/pnProjX[cIdx]):0;
    for (cIdx=R/10;(pnProjL_s[cIdx+1]<=pnProjL_s[cIdx])&&(cIdx<R);cIdx++);
    for (;cIdx>=0;pnProjL_s[cIdx--]=1);
    for (cIdx=R*3/4;(pnProjL_s[cIdx-1]<=pnProjL_s[cIdx])&&(cIdx);cIdx--);
    for (;cIdx<R;pnProjL_s[cIdx++]=1);
memcpy(pnProjL,pnProjL_s,R*4);
    IPL_HIST_Blur(pnProjX_s,pnProjB_s,R,BLWIN_I2P);
    IPL_HIST_Blur(pnProjX,pnProjB,R,BLWIN_I2P);
    for (cIdx=0;cIdx<R;cIdx++)
      pnProjB_s[cIdx] = (pnProjX[cIdx])?(1000*pnProjX_s[cIdx]/pnProjX[cIdx]):0;
    for (cIdx=R/10;(pnProjB_s[cIdx+1]<=pnProjB_s[cIdx])&&(cIdx<R);cIdx++);
    for (;cIdx>=0;pnProjB_s[cIdx--]=1);
    for (cIdx=R*3/4;(pnProjB_s[cIdx-1]<=pnProjB_s[cIdx])&&(cIdx);cIdx--);
    for (;cIdx<R;pnProjB_s[cIdx++]=1);
memcpy(pnProjB,pnProjB_s,R*4);
    IPL_HIST_Blur(pnProjX_s,pnProjT_s,R,BLWIN_I2P);
    IPL_HIST_Blur(pnProjX,pnProjT,R,BLWIN_I2P);
    for (cIdx=0;cIdx<R;cIdx++)
      pnProjT_s[cIdx] = (pnProjX[cIdx])?(1000*pnProjX_s[cIdx]/pnProjX[cIdx]):0;
    for (cIdx=R/10;(pnProjT_s[cIdx+1]<=pnProjT_s[cIdx])&&(cIdx<R);cIdx++);
    for (;cIdx>=0;pnProjT_s[cIdx--]=1);
    for (cIdx=R*3/4;(pnProjT_s[cIdx-1]<=pnProjT_s[cIdx])&&(cIdx);cIdx--);
    for (;cIdx<R;pnProjT_s[cIdx++]=1);
memcpy(pnProjT,pnProjT_s,R*4);
    // free & exit
    if (dbgim)
      free(dbgim);
    return res;
}
*/
