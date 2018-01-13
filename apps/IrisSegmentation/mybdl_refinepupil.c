#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "header.h"
#include "imio.h"
#include "ipl.h"
#include "IrisSegmentation.h"

#define Cder 2
#define MAXCOST 214748364

typedef struct
{
	int cost;
	int* path;
} SCirCand;

typedef struct
{
	int cost, root, rad;
} SPathCand;

static int sort_PC(const void* e1, const void* e2)
{
    return ((SPathCand*)e1)->cost-((SPathCand*)e2)->cost;
}


/*int CircBypass(int start, unsigned char* priceMap, int* radPath, int pathCost, int H, int W)
{
	//FILE *fout = fopen("../data/res/0612/EnhancementMaps.txt","w+");
	SPathCand* PC;
	int i,n, sPnt, costUpd, c1, c2, c3, res, maxprice, nC, optRoot, optRad, optCost;
	int *rootMap, *costMap, *parMap;
	rootMap = (int*)malloc(H*W*sizeof(int));
	costMap = (int*)malloc(H*W*sizeof(int));
	parMap = (int*)malloc(H*W*sizeof(int));

	maxprice = 2*H*W;
	memset(costMap,maxprice,H*W*sizeof(int));
	memset(rootMap,-1,H*W*sizeof(int));
	memset(parMap,-1,H*W*sizeof(int));
	
	res = 0;
	for (i = 0; i < H; i++)
	{
		rootMap[i*W] = i;
		parMap[i*W] = i;
		if (i < start + 3 && i > start - 3)
			costMap[i*W] = priceMap[i*W];
		else
		{
			costMap[i*W] = maxprice;//maximum price for passing all the graph
			res++;
		}
	}
	if (res == H)
	{
		free(rootMap);
		free(costMap);
		free(parMap);
		return -1;
	}

	for (n = 1; n < W; n++)
	{
		c1 = costMap[n-1] + (int)priceMap[n-1];
		c2 = costMap[W + n-1] + (int)priceMap[W + n-1] + Cder;//Изменить потом величину штрафа за величину производной
		if (c1 < c2)
		{
			if (c1 >= 0)
			{
				costMap[n] = c1;
				rootMap[n] = rootMap[n-1];
				parMap[n] = 0;
			}
			else
			{
				costMap[n] = c2;
				rootMap[n] = rootMap[n-1 + W];
				parMap[n] = 1;
			}
		}
		else
		{
			if (c2 >= 0)
			{
				costMap[n] = c2;
				parMap[n] = 1;
				rootMap[n] = rootMap[W + n-1];
			}
			else
			{
				costMap[n] = c1;
				rootMap[n] = rootMap[n-1];
				parMap[n] = 0;
			}
		}

		for (i = 1; i < H-1; i++)
		{
			c1 = costMap[i*W + n-1] + (int)priceMap[i*W + n-1];
			c2 = costMap[(i+1)*W + n-1] + (int)priceMap[(i+1)*W + n-1] + Cder;
			c3 = costMap[(i-1)*W + n-1] + (int)priceMap[(i-1)*W + n-1] + Cder;

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
						sPnt = i-1;
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
						sPnt = i+1;
					}
				}
				else
				{
					if (c3 >= 0)
					{
						costUpd = c3;
						sPnt = i-1;
					}
				}
			}

			costMap[i*W + n] = costUpd;
			rootMap[i*W + n] = rootMap[sPnt * W + n-1];
			parMap[i*W + n] = sPnt;
		}

		c1 = costMap[(H-1) * W + n-1] + (int)priceMap[(H-1) * W + n-1];
		c2 = costMap[(H-2) * W + n-1] + (int)priceMap[(H-2) * W + n-1] + Cder;//Изменить потом величину штрафа за величину производной
		if (c1 < c2)
		{
			if (c1 >= 0)
			{
				costMap[(H-1)*W + n] = c1;
				parMap[(H-1)*W + n] = H-1;
				rootMap[(H-1)*W + n] = rootMap[(H-1)*W + n-1];
			}
			else
			{
				costMap[(H-1)*W + n] = c2;
				parMap[(H-1)*W + n] = H-2;
				rootMap[(H-1)*W + n] = rootMap[(H-2)*W + n-1];
			}
		}
		else
		{
			costMap[(H-1)*W + n] = c2;
			parMap[(H-1)*W + n] = H-2;
			rootMap[(H-1)*W + n] = rootMap[(H-2)*W + n-1];
		}
	}


	PC = (SPathCand*)malloc(H * sizeof(SPathCand));
	nC = 0;
	if (start < 2)
		sPnt = 0;
	else
		sPnt = start - 2;
	for (i = sPnt; i < start + 3 && i < H; i++)
	{
		if (rootMap[i * W + W-1] > start + 2 || rootMap[i * W + W-1] < start - 2)
			costMap[i * W + W-1] = maxprice;
		else
		{
			PC[nC].cost = costMap[i*W + W-1];
			PC[nC].rad = i;
			PC[nC].root = rootMap[i*W + W-1];
			nC++;
		}
	}


	if (nC == 0)
		optCost = -1;
	else
	{
		qsort(PC,nC,sizeof(PC[0]),sort_PC);
		optRad = PC[0].rad;
		optRoot = PC[0].root;
		optCost = PC[0].cost;
		if (optCost > 0)
		{
			radPath[W] = (optRoot + optRad)/2;
			radPath[W-1] = optRad;
			for (i = W-2; i >=0; i--)
				radPath[i] = parMap[i + radPath[i+1]*W];
		}
		
	}
	free(rootMap);
	free(costMap);
	free(parMap);
	free(PC);
	return optCost;
	
}

int pupilCircSP(int* radPath, unsigned char *Map, int H, int W, int pRad)
{
	int res;
	int r, i,j, cC, tmp, argmin, min = H*W;
	unsigned char t;
	int *cand;
	SCirCand *CSPath;
	cand = (int*)malloc(H * sizeof(int));
	CSPath = (SCirCand*)malloc(H * sizeof(SCirCand));
	cC = 0;
	for (i = 0; i < H; i++)
	{
		t = Map[i*W];
		if (t == 0)
		{
			cand[cC] = i;
			CSPath[cC].path = (int*)malloc(361 * sizeof(int));
			CSPath[cC].path[0] = i;
			CSPath[cC].cost = Map[i * W];
			cC++;
		}
	}
	printf("CC=%d\n", cC);
	argmin = -1;
	res = 0;
	if (cC > 0)
	{
		for (i = 0; i < cC; i++)
		{
			//printf("%d ", cand[i]);
			tmp = CircBypass(cand[i], Map, CSPath[i].path, CSPath[i].cost,H,W);
			if (tmp > 0 && tmp < min)
			{
				min = tmp;
				argmin = i;
			}
			else
			{
				if (argmin == -1)
				{
					res = -1;
					continue;
				}
				else
				{
					res = 0;
				}
			}
			CSPath[i].cost = tmp;
		}
		//printf("\n");
		if (res != -1)
		{
			for (j = 0; j < W+1; j++)
				radPath[j] = CSPath[argmin].path[j];
		}
		else
		{
			for(j = 0; j < cC; j++)
				free(CSPath[j].path);
			free(CSPath);
			free(cand);
			return -1;
		}
		
	}
	else
	{
		CSPath[0].path = (int*)malloc(361 * sizeof(int));
		CSPath[0].path[0] = pRad;
		tmp = CircBypass(pRad, Map, CSPath[0].path, CSPath[0].cost,H,W);

		for (j = 0; j < W+1; j++)
			radPath[j] = CSPath[0].path[j];
	}
	for(j = 0; j < cC; j++)
		free(CSPath[j].path);
	free(CSPath);
	free(cand);
	return 0;
}


int IBD_RefinePupil(int *xpList,int *ypList, int* radmean, CInfo *CP, int MinRad, int MaxRad, char *name, const unsigned char* img, int H, int W, int flags)
{
	unsigned char *imgEdge, *imgCrop, *imgGrad;
	void *canny_buf = NULL;
	int Wseg, Hseg, i, rad, angle, xbeg,ybeg,xend,yend, Hc, Wc, x, y, xc,yc, med, canny_sz,res;
	int gx,gy, ax,ay, *radList, max;
	int *dest, *imdx, *imdy;
	double t, tlow, thigh, sigma, val;
	int16 *pgx_calc,*pgy_calc;
	//FILE *fout = fopen("../data/res/edgeMapVals.txt", "w");

	//Cropping the pupil for better performance
	xbeg = CP->xc - CP->r - W/10;
	if (xbeg < 0)	xbeg = 1;
	ybeg = CP->yc - CP->r - H/10;
	if (ybeg < 0)	ybeg = 1;
	xend = CP->xc + CP->r + W/10;
	if (xend > W-1)	xend = W-1;
	yend = CP->yc + CP->r + H/10;
	if (yend >= H)	yend = H;
	if (yend < 0 || xbeg >= W || ybeg >= H || xend < 0)	return -1;
	Hc = yend-ybeg + 1;
	Wc = xend - xbeg + 1;
	xc = CP->xc - xbeg;
	yc = CP->yc - ybeg;

	imgCrop = (unsigned char*)malloc(Hc * Wc * sizeof(unsigned char));
	imgGrad = (unsigned char*)malloc(Hc * Wc * sizeof(unsigned char));
	imgEdge = (unsigned char*)malloc(Hc * Wc * sizeof(unsigned char));
	dest = (int*)malloc(Hc * Wc * sizeof(int));
	imdx = (int*)malloc(Hc * Wc * sizeof(int));
	imdy = (int*)malloc(Hc * Wc * sizeof(int));
	
	memset(dest,0,Hc*Wc*sizeof(int));
	memset(imgCrop,0,Hc*Wc);
	memset(imgGrad,0,Hc*Wc);
	canny_buf = (void*)malloc(20 * H * W*sizeof(unsigned char));

	for (y = ybeg; y <= yend; y++)
	{
		for (x = xbeg; x <= xend; x++)
		{
			imgCrop[(y-ybeg)* Wc + (x - xbeg)] = img[y * W + x];
		}
	}

	med = IPL_HIST_mdn_pixE(imgCrop, Hc*Wc);
	//printf("med_crop=%d\n", med);
	t = (double)med/255;
	//printf("t = %f\n", t);
	tlow =  0.8 * t;
	thigh = 1.33 * t;
	//tlow = .55f;
	//thigh = .7f;
	sigma = 2.f + 3*t;
	canny_sz = 20 * H * W;//1049665; //Костыль

	//FLT_Sobel3x3(dest, imgEdge, imgCrop, name, Hc, Wc, mask, imdx, imdy,(FLT_SOBEL_SELECTION));


	if ((res = IPL_FILT_Canny(imgEdge,imgCrop,&pgx_calc,&pgy_calc,Wc,Hc,sigma,tlow,thigh,
				canny_buf,&canny_sz,name))!=ERR_OK)					//original: 4.f .7f .8f
	{
	  printf("Error %d: Canny, refinement stage fault\n", res);	
	  free(imgCrop);
	  free(imgGrad);
	  free(imgEdge);
	  free(canny_buf);
	  return -1;		
	}


	max = 0;
	for (y = 1; y < Hc-1; y++)
	{
		for (x = 1; x < Wc-1; x++)
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
			t = 180. * acos(val)/PI;
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

	Wseg = 360;
	Hseg = MaxRad - MinRad + 1;
	radList = (int*)malloc((Wseg+1) * sizeof(int));
	printf("improving\n");
	
	imgPolar = (unsigned char*)malloc(Hseg*Wseg*sizeof(unsigned char));
	memset(imgPolar,0,Hseg*Wseg*sizeof(unsigned char));
	

	for (rad = MinRad; rad <= MaxRad; rad++)
	{
		for (angle = 0; angle < Wseg; angle++)
		{
			i = (int)((double)rad * sin((double)(angle)/180 * PI)) * Wc  + (int)((double)rad * cos((double)(angle)/180 * PI));
			imgPolar[(rad-MinRad)*Wseg + angle] = 255 - imgEdge[Wc * (CP->yc - ybeg) + (CP->xc - xbeg) + i];

		}
	}	
	for (rad = 0; rad < Hseg; rad++)
	{
		for (angle = 0; angle < Wseg; angle++)
			if (imgPolar[rad * Wseg + angle] > 0)
				imgPolar[rad * Wseg + angle] = (unsigned char)1;
	}
	if ((res = pupilCircSP(radList, imgPolar, Hseg, Wseg, CP->r - MinRad)) != 0)
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
	*radmean = 0;
	for (i = 0; i < Wseg; i++)
	{
		radList[i] += MinRad;
		*radmean += radList[i];
		//printf("%d ", radList[i]);
		xpList[i] = CP->xc + (int)((double)radList[i]*cos((double)(i)/180 * PI));
		ypList[i] = CP->yc + (int)((double)radList[i]*sin((double)(i)/180 * PI));
	}
	*radmean /= Wseg;

	free(radList);
	free(canny_buf);
	free(imgCrop);
	free(imgGrad);
	free(imgEdge);
	free(dest);
	free(imdx);
	free(imdy);
	return 0;
}*/