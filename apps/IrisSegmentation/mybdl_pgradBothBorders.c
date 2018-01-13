#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "header.h"
#include "imio.h"
#define ANGDELT 17
#define RADDELT 2
#define MINRAD 20

typedef struct
{
  int32 dir;
  int16 x;
  int16 y;
  int16 gx;
  int16 gy;
} SEdgePnt;


typedef struct
{
	int16 xc;
	int16 yc;
	int16 r;
	int16 xm;
	int16 ym;
} SCircleCnd;

static int sort_edgepnt(const void* e1, const void* e2)
{
    return ((SEdgePnt*)e1)->dir-((SEdgePnt*)e2)->dir;
}

static int sort_img8U(const void* e1, const void* e2)
{
	return *((unsigned char*)e1) - *((unsigned char*)e2);
}


int ibd_graparBothBorders(CInfo* CI, CInfo* CP, char* name, const unsigned char* img, unsigned char *imgBl, int H, int W, int dPhi, int flags)
{	
	int nGp,cGp,gval,val,scalar, x,y,x0,y0,cPt,nBegPairIdx,nEndPairIdx,tint,cPar,D,Dt,rmn,rmx, nC,cC, vtype = 0;
	int dx,dy,dr1,dr2;
	double t, xt,yt,r1,r2,dx1,dx2,dy1,dy2,xc,yc,dr ,sigma,tlow,thigh, vrmax;
	int xmax=-1,ymax=-1,vmax,rmax=-1, maxRad;
	uint8 med, max;
	FILE* rhout;

	SEdgePnt* pnEdgePnts;
	SCircleCnd* pnCircleCnd;
	int *pnAccum,*pnAccBl;
	uint8* pucEdge=NULL;
	uint8 *imgsort, *imgMorph;
	void* canny_buf=NULL;
	double* radHist, *radHistSmooth;
	int16 *pgx_calc,*pgy_calc;
	
	int canny_sz;
	int sz;
	int res = 0;
    sz = 0;
    sz += H*W*sizeof(pnEdgePnts[0])+  // edge array
          H*W*sizeof(pnAccum[0])+     // accum
          H*W*sizeof(pnAccBl[0])+     // bluraccum
		  H*W*sizeof(pnCircleCnd[0]); // circle params 
	

    pnEdgePnts = (SEdgePnt*)malloc(H*W*sizeof(pnEdgePnts[0]));
    pnAccum = (int*)malloc(H*W*sizeof(pnAccum[0]));
	pnAccBl = (int*)malloc(H*W*sizeof(pnAccBl[0]));
	pnCircleCnd = (SCircleCnd*)malloc(2*H*W*sizeof(pnCircleCnd[0]));
	pucEdge = (uint8*)malloc(H*W*sizeof(pucEdge[0]));
	canny_buf = (void*)malloc(4*sz);

	//median search for canny tresholding
	imgsort = (unsigned char*)malloc(H*W*sizeof(imgsort[0]));
	memcpy(imgsort,img,sizeof(img[0])*H*W);
	qsort(imgsort,H*W,sizeof(uint8),sort_img8U);
	if (H*W % 2 == 1)
		med = imgsort[H*W/2 + 1];
	else
		med = (imgsort[H*W/2] + imgsort[H*W/2])/2;
	max = imgsort[H*W-1];
	free(imgsort);

	//Canny edge detection and intensity gradients calculating
	//automated canny tresholds calc
	t = (double)med/(int)max;
	printf("t = %f, max = %u\n", t,max);
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
	if ((res = IPL_FILT_Canny(pucEdge,imgBl,&pgx_calc,&pgy_calc,W,H,
				sigma, //5.f,
				tlow,
				thigh,
				canny_buf,&canny_sz,name))!=ERR_OK)					//original: 4.f .7f .8f
	{
	  printf("Error %d: Canny\n", res);	
	  //return -1;		
	}

	if (flags&BDL_PGRAD_SAVECANNY)
		SaveBmp8(name,"_canny.bmp", W, H, pucEdge, 4);


	// Collect gradient point to array
    nGp = 0;
	for (y = 10; y < H-10; y++)
	{
		for (x = 10; x < W-10; x++)
		{

			if (pucEdge[cPt = y*W + x])
			{
				pnEdgePnts[nGp].dir = (int)(.5+180.*atan2((double)pgy_calc[cPt],(double)pgx_calc[cPt])/PI);
				pnEdgePnts[nGp].x = (int16)x;
				pnEdgePnts[nGp].y = (int16)y;
				pnEdgePnts[nGp].gx = pgx_calc[cPt];
				pnEdgePnts[nGp].gy = pgy_calc[cPt];
				nGp++;
			}
		}
	}
	
	
	//sorting the gradient points array
	qsort(pnEdgePnts,nGp,sizeof(pnEdgePnts[0]),sort_edgepnt);
	
	//duplicate the array
	memcpy(pnEdgePnts+nGp,pnEdgePnts,sizeof(pnEdgePnts[0])*nGp);
	
	//clean the accumulator
    for (cGp=0;cGp<nGp;cGp++)
		pnEdgePnts[cGp+nGp].dir = pnEdgePnts[cGp].dir+360;
		
    // clean separated accumulators for center and radius
    memset(pnAccum, 0, W*H*sizeof(pnAccum[0]));

	// prime counters
    tint = pnEdgePnts[0].dir+dPhi-ANGDELT;
    for (nBegPairIdx=0;((nBegPairIdx<nGp)&&
          (pnEdgePnts[nBegPairIdx].dir<tint));nBegPairIdx++);

    tint = pnEdgePnts[0].dir+dPhi+ANGDELT;
    for (nEndPairIdx=nBegPairIdx;((nEndPairIdx<nGp)&&
          (pnEdgePnts[nEndPairIdx].dir<tint));nEndPairIdx++);
	
	if (flags&BDL_PGRAD_LINEVOTE)
		vtype = 1;

	nC = 0;
    // main processing
	for (cGp=0;cGp<nGp;cGp++)
	{
	  for (cPar=nBegPairIdx + cGp; cPar < nEndPairIdx + cGp; cPar++)
	  {
		D = pnEdgePnts[cPar].gx * pnEdgePnts[cGp].gy - pnEdgePnts[cPar].gy * pnEdgePnts[cGp].gx;
		Dt = pnEdgePnts[cPar].gx * (pnEdgePnts[cPar].y-pnEdgePnts[cGp].y) - pnEdgePnts[cPar].gy * (pnEdgePnts[cPar].x-pnEdgePnts[cGp].x);
		
		if (D != 0)
		{  //continue;
			t = (double)Dt / D;
			dx1 = (double)pnEdgePnts[cGp].gx * t;
			dy1 = (double)pnEdgePnts[cGp].gy * t;
			xt = (double)pnEdgePnts[cGp].x + dx1;
			yt = (double)pnEdgePnts[cGp].y + dy1;
			if (xt < 1 || yt < 1 || xt >= W-1 || yt >= H-1)
				continue;
			dx2 = (double)pnEdgePnts[cPar].x - xt;
			dy2 = (double)pnEdgePnts[cPar].y - yt;
			r1 = sqrt( dx1*dx1 + dy1*dy1 );
			r2 = sqrt( dx2*dx2 + dy2*dy2 );
			if (abs(r1 - r2) >= RADDELT || max(r1,r2) < MINRAD)
				continue;
			//Ds = pnEdgePnts[cGp].gx * (pnEdgePnts[cPar].y-pnEdgePnts[cGp].y) - pnEdgePnts[cGp].gy * (pnEdgePnts[cPar].x-pnEdgePnts[cGp].x);
			//s = (double)Ds / (double)D;
			pnCircleCnd[nC].xc = (int)xt;
			pnCircleCnd[nC].yc = (int)yt;
			pnCircleCnd[nC].r = (r1+r2)/2;
			pnCircleCnd[nC].xm = 0;
			pnCircleCnd[nC].ym = 0;
			nC++;
			switch(vtype) {
				case 0:
					x = (int)xt;
					y = (int)yt;
					pnAccum[y*W+x]++;
					break;
				case 1:
					dx1 = 0.5 * (double)(pnEdgePnts[cGp].gx + pnEdgePnts[cPar].gx);
					dy1 = 0.5 * (double)(pnEdgePnts[cGp].gy + pnEdgePnts[cPar].gy);
					x = (int)(xt + dx1/sqrt(dx1*dx1 + dy1*dy1)*5);
					y = (int)(yt + dy1/sqrt(dx1*dx1 + dy1*dy1)*5);
					x0 = (int)(xt + dx1/sqrt(dx1*dx1 + dy1*dy1)*2);
					y0 = (int)(yt + dy1/sqrt(dx1*dx1 + dy1*dy1)*2);
					DRAW_Line(pnAccum, H, W, x, y, x0,y0);
					x = (int)(xt - dx1/sqrt(dx1*dx1 + dy1*dy1)*5);
					y = (int)(yt - dy1/sqrt(dx1*dx1 + dy1*dy1)*5);
					x0 = (int)(xt - dx1/sqrt(dx1*dx1 + dy1*dy1)*2);
					y0 = (int)(yt - dy1/sqrt(dx1*dx1 + dy1*dy1)*2);
					DRAW_Line(pnAccum, H, W, x, y, x0,y0);
					break;
			}
		}
	  }
	}

	maxRad = min(H,W)/2;
	// blur accumulator
    IPL_FILT_HaussBlur3x3_Int(pnAccBl,pnAccum,W,H,0);
	memcpy(pnAccum, pnAccBl, W*H*sizeof(int));
	IPL_FILT_HaussBlur3x3_Int(pnAccBl,pnAccum,W,H,0);

	// find maximum
    vmax = xmax = ymax = 0;
    for (y=1;y<H-1;y++)
      for (x=1;x<W-1;x++)
        if (vmax < pnAccBl[y*W+x])
          vmax = pnAccBl[(ymax = y)*W+(xmax = x)];
	CI->xc = xmax+1;
	CI->yc = ymax+1;
	if (vmax > 0)
	{
		unsigned char* pnAccOut = (unsigned char*)malloc(H*W*sizeof(unsigned char));
		// Cleaning borders for accumulator
		for(x = 0; x < max(W,H); x++)
		{
			pnAccOut[x] = 0;
			pnAccOut[(H-1)*W + x] = 0;
			if (x < min(W,H))
			{
				pnAccOut[x*(W-1)] = 0;
				pnAccOut[x*(W-1) + H-1] = 0;
			}
		}
		// Intensity normalization
		for(y = 1; y < H-1; y++)
			for(x = 1; x < W-1; x++)
				pnAccOut[y*W + x] = (unsigned char)((double)pnAccBl[y*W + x] / vmax * 255);

		if (flags&BDL_PGRAD_SAVEACC)
		{
			SaveBmp8(name, "_acc1.bmp", W, H, pnAccOut, 4);
		}
		free(pnAccOut);
	}

	vrmax = 0;
	vmax = rmax = 0;
	if(flags&BDL_PGRAD_CALCRADIUS)
	{
		rmx = max(H,W)/2;
		radHist = (double*)malloc((H+W)*sizeof(double));
		radHistSmooth = (double*)malloc((H+W)*sizeof(double));
		for(x = 0; x < 400; x++) radHist[x] = 0.0;
		
		for(cGp=0;cGp<nGp;cGp++)
		{
			dx1 = (double)(pnEdgePnts[cGp].x - xmax);
			dy1 = (double)(pnEdgePnts[cGp].y - ymax);
			xt = (double)pnEdgePnts[cGp].gx;
			yt = (double)pnEdgePnts[cGp].gy;
			t = (dx1*pnEdgePnts[cGp].gx + dy1 * pnEdgePnts[cGp].gy)/ sqrt(xt * xt + yt * yt);
			r1 = sqrt( dx1*dx1 + dy1*dy1 );
			t = abs(t) / r1;
			if(r1 >= MINRAD && r1 < rmx && t > 0.96)
				radHist[(int)(r1+.5)]++;
		}

		// Histogram normalization
		for (x = 1; x < rmx; x++) radHist[x] = radHist[x]/x; 
		//Histogram Haussian filtering (twice for accuracy)
		IPL_HIST_Blur_double(radHistSmooth, (const double*) radHist, rmx, 5);
		for (x = 1; x < rmx; x++) radHist[x] = radHistSmooth[x];
		IPL_HIST_Blur_double(radHistSmooth, (const double*) radHist, rmx, 5);


		/*rhout = fopen("../data/res/radHist1.txt", "wt");
		for(x = 1; x < rmx-2; x++){
			fprintf(rhout, "%1.4f  ", radHistSmooth[x]);
		}
		fprintf(rhout, "\n");
		fclose(rhout);*/
		
		for (x=MINRAD; x<rmx; x++)
			if (vrmax < radHistSmooth[x])
				vrmax = radHistSmooth[rmax = x];
		CI->r = rmax;
	}


	//Second border search

	//clean accumulator
	memset(pnAccum, 0, W*H*sizeof(pnAccum[0]));

	dr2 = CI->r / 5;
	dr2 = dr2*dr2;
	//CI->yc = H - CI->yc;
	for (cC=0; cC < nC; cC++)
	{
		dx = pnCircleCnd[cC].xc - CI->xc;
		dy = pnCircleCnd[cC].yc - CI->yc;
		dr1 = dx*dx + dy*dy;
		if (dr1 <= dr2 && (pnCircleCnd[cC].r < 3*CI->r / 4 ||  pnCircleCnd[cC].r > 4*CI->r / 3))
		{
			pnAccum[pnCircleCnd[cC].yc*W+pnCircleCnd[cC].xc]++;
			//radHist[pnCircleCnd[cC].r] *= pnCircleCnd[cC].r;
			//radHistSmooth[pnCircleCnd[cC].r]++;
			//radHist[pnCircleCnd[cC].r] /= pnCircleCnd[cC].r;
			//printf("%d \n", pnCircleCnd[cC].r);
		}
	}
	// blur accumulator
    IPL_FILT_HaussBlur3x3_Int(pnAccBl,pnAccum,W,H,0);
	memcpy(pnAccum,pnAccBl, H*W*sizeof(int));
	IPL_FILT_HaussBlur3x3_Int(pnAccBl,pnAccum,W,H,0);

	// find maximum
    vmax = xmax = ymax = 0;
    for (y=10;y<H-10;y++)
	{
		for (x=10;x<W-10;x++)
	  {
		  if (vmax < pnAccBl[y*W+x])
			vmax = pnAccBl[(ymax = y)*W+(xmax = x)];
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
	
	if(flags&BDL_PGRAD_CALCRADIUS)
	{

		for (dr1 = 0; dr1 <= CI->r / 7; dr1++)
			radHistSmooth[dr1] = 0.0;
		for (dr1 = 0; dr1 <= MINRAD; dr1++)
			radHistSmooth[dr1] = 0.0;
		for (dr1 = 3 * CI->r/4; dr1 < 4 * CI->r /3; dr1++)
			radHistSmooth[dr1] = 0.0;
	
		rmx = max(H,W)/2;
		
		for (dr1 = 0; dr1 < rmx; dr1++)
			radHist[dr1] = radHistSmooth[dr1];
		IPL_HIST_Blur_double(radHistSmooth, (const double*) radHist, rmx, 3);
		//'Haussian' filtering for radius search
		/*memcpy(radHist, radHistSmooth, sizeof(radHist[0]));
		for(x = 3; x < 398; x++)
			radHistSmooth[x] = 0.0938 * radHist[x-4] + 0.1061 * radHist[x-3] + 0.1158 * radHist[x-2] + 0.1221 * radHist[x-1] +
				0.1243 * radHist[x] + 0.1221 * radHist[x+1] + 0.1158 * radHist[x+2] + 0.1061 * radHist[x+3] + 0.0938 * radHist[x+4];
		memcpy(radHist, radHistSmooth, sizeof(radHist[0]));
		*/
		/*rhout = fopen("../data/res/radHist2.txt", "wt");
		for(x = 3; x < 398; x++){
			fprintf(rhout, "%1.2f  ", radHistSmooth[x]);
		}
		
		fprintf(rhout, "\n");
		fclose(rhout);
		*/
		vrmax = 0.0;
		for (x=rmx-1; x>MINRAD; x--)
			if (vrmax < radHistSmooth[x])
				vrmax = radHistSmooth[rmax = x];
		CP->r = rmax;
		free(radHist);
		free(radHistSmooth);
	}

	if (CP->r > CI->r)
	{
		val = CP->r;
		CP->r = CI->r;
		CI->r = val;
		val = CP->xc;
		CP->xc = CI->xc;
		CI->xc = val;
		val = CP->yc;
		CP->yc = CI->yc;
		CI->yc = val;
	}
	free(pnEdgePnts);
	free(pnAccum);
	free(pnAccBl);
	free(pnCircleCnd);
	free(pucEdge);
	free(canny_buf);
	return 0;
}