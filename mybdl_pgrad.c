#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "header.h"
#include "imio.h"
#define ANGDELT 17
#define RADDELT 3
#define MINRAD 20


typedef struct
{
  int32 dir;
  int16 x;
  int16 y;
  int16 gx;
  int16 gy;
} SEdgePnt;





static int sort_edgepnt(const void* e1, const void* e2)
{
    return ((SEdgePnt*)e1)->dir-((SEdgePnt*)e2)->dir;
}

int ibd_grapar(CInfo* CI, char* name, const unsigned char* img, int H, int W, int dPhi, SBoundPnt* pnBoundary, CInfo* CP, int flags)
{	
	char* cname;
	int nGp,cGp,nBp,gval, val, scalar, x,y,x0,y0,cPt,nBegPairIdx,nEndPairIdx,tint,cPar,D,Dt,rmn,rmx, vtype = 0;
	double t, xt,yt,r1,r2,dx1,dx2,dy1,dy2,xc,yc,dr,sigma,tlow,thigh;
	int xmax=-1,ymax=-1,vmax,rmax=-1, Tg=430;
	SEdgePnt* pnEdgePnts;

	int *pnAccum,*pnAccBl, *pnRadHist;
	uint8* pucEdge=NULL;
	void* canny_buf=NULL;
	int canny_sz;
	int sz;
	int res = 0;
	//char nam[FILENAME_MAX];
	int16 *pgx_calc,*pgy_calc;

    unsigned char *img_bl8 = (unsigned char*)malloc(H*W*sizeof(unsigned char));
	IPL_FILT_HaussBlur5x5(img_bl8,img,W,H);
	IPL_FILT_HaussBlur5x5(img_bl8,img_bl8,W,H);
	IPL_FILT_HaussBlur5x5(img_bl8,img_bl8,W,H);
	IPL_FILT_HaussBlur5x5(img_bl8,img_bl8,W,H);
	IPL_FILT_HaussBlur5x5(img_bl8,img_bl8,W,H);
	


    sz = 0;
    sz += H*W*sizeof(pnEdgePnts[0])+  // edge array
          H*W*sizeof(pnAccum[0])+     // accum
          H*W*sizeof(pnAccBl[0])+     // bluraccum
          (H+W)*sizeof(pnRadHist[0]); // radhist 

	void* pvBuf = (void*)malloc(3 * sz);//more than we need

    pnEdgePnts = (SEdgePnt*)pvBuf;
    pnAccum = (int*)(pnEdgePnts + H*W);
    pnAccBl = pnAccum + H*W;
    pnRadHist = pnAccBl + H*W;
	
	pucEdge = (uint8*)(pnRadHist+H+W);
	canny_buf = (void*)(pucEdge+H*W);
	if (CP == NULL)
	{
		tlow = .6f;
		thigh = .8f;
		sigma = 5.f;
	}
	else
	{
		tlow = .4f;
		thigh = .6f;
		sigma = 7.f;
	}

	canny_sz = 5222400;
	if ((res = IPL_FILT_Canny(pucEdge,img_bl8,&pgx_calc,&pgy_calc,W,H,
				sigma, //5.f,
				tlow,
				thigh,
				canny_buf,&canny_sz,name))!=ERR_OK)					//original: 4.f .7f .8f
	{
	  printf("Error %d: Canny\n", res);	
	  return -1;		
	}
	free(img_bl8);

	if (CP != NULL)
	{
		xc = (double)CP->xc;
		yc = (double)CP->yc;
		rmn = (int)(0.75*0.75 * CP->r*CP->r);
		rmx = (int)(16 * (double)(CP->r * CP->r)/9);
		dr = 10;
		for (y = 1; y < H-1; y++)
			for (x = 1; x < W-1; x++)
				if (pucEdge[cPt = y*W + x])
				{
					int dr = (x - xc) * (x - xc) + (y - yc) * (y - yc);
					if (dr > rmn & dr < rmx)
						pucEdge[cPt] = 0;
				}
	}
	if (flags&BDL_PGRAD_SAVECANNY)
	{
		cname = (char*)malloc(256 * sizeof(char));
		strcpy(cname, name);
		if (CP == NULL)
			strcat(cname, "_canny.bmp");
		else
			strcat(cname, "_canny2.bmp");
		CreateBmp8(cname, W, H, pucEdge, 4);
		free(cname);
	}
	//FILE* out = fopen("pgy_calc.txt", "w");
	
	// collect gradient point to array
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
			//fprintf(out, "%d ", pucEdge[cPt]);
		}
		//fprintf(out, "\n");
	}
	//fclose(out);
	
	
	//sorting the gradient points array
	qsort(pnEdgePnts,nGp,sizeof(pnEdgePnts[0]),sort_edgepnt);
	//duplicate the array
	memcpy(pnEdgePnts+nGp,pnEdgePnts,sizeof(pnEdgePnts[0])*nGp);
	//clean the accumulator
    for (cGp=0;cGp<nGp;cGp++)
      pnEdgePnts[cGp+nGp].dir = pnEdgePnts[cGp].dir+360;
    // clean accumulator
    memset(pnAccum, 0, W*H*sizeof(pnAccum[0]));
	memset(pnRadHist,0,(W+H)*sizeof(pnRadHist[0]));

	// prime counters
    tint = pnEdgePnts[0].dir+dPhi-ANGDELT;
    for (nBegPairIdx=0;((nBegPairIdx<nGp)&&
          (pnEdgePnts[nBegPairIdx].dir<tint));nBegPairIdx++);

    tint = pnEdgePnts[0].dir+dPhi+ANGDELT;
    for (nEndPairIdx=nBegPairIdx;((nEndPairIdx<nGp)&&
          (pnEdgePnts[nEndPairIdx].dir<tint));nEndPairIdx++);
	

	if (flags&BDL_PGRAD_LINEVOTE)
		vtype = 1;

    // main processing
	if (CP == NULL)
	{
		for (cGp=0;cGp<nGp-nEndPairIdx;cGp++)
		{
		  for (cPar=nBegPairIdx+cGp;cPar<nEndPairIdx+cGp;cPar++)
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
	}
	else
	{
		for (cGp=0;cGp<nGp;cGp++)
		{
		  for (cPar=nBegPairIdx+cGp;cPar<nEndPairIdx+cGp;cPar++)
		  {
			D = pnEdgePnts[cPar].gx * pnEdgePnts[cGp].gy - pnEdgePnts[cPar].gy * pnEdgePnts[cGp].gx;
			Dt = pnEdgePnts[cPar].gx * (pnEdgePnts[cPar].y-pnEdgePnts[cGp].y) - pnEdgePnts[cPar].gy * (pnEdgePnts[cPar].x-pnEdgePnts[cGp].x);
			
			if (D == 0)
			  continue;
			t = (double)Dt / (double)D;
			//printf("%f\n", t);
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
			if (xt < xc + dr & xt > xc - dr & yt > yc - dr & yt < yc + dr)
			{
				x = (int)xt;
				y = (int)yt;
				pnAccum[y*W+x]++;
			}
			//pnRadHist[(int)(r1+.5)]++;
		  }
		}
		

	}
	// blur accumulator
    pnAccBl = (int*)malloc(W*H*sizeof(pnAccBl[0]));
    IPL_FILT_HaussBlur3x3_Int(pnAccBl,pnAccum,W,H,0);
	IPL_FILT_HaussBlur3x3_Int(pnAccBl,pnAccBl,W,H,0);

	// find maximum
    vmax = xmax = ymax = 0;
    for (y=1;y<H-1;y++)
      for (x=1;x<W-1;x++)
        if (vmax < pnAccBl[y*W+x])
          vmax = pnAccBl[(ymax = y)*W+(xmax = x)];
	
	if (CP != NULL)
	{
		if( (xmax - xc)*(xmax - xc) + (ymax - yc) * (ymax - yc) >= (0.15*0.15*CP->r*CP->r))
		{
			xmax = xc;
			ymax = yc;
		}

	}
	if (vmax != 0)
	{
		unsigned char* pnAccOut = (unsigned char*)malloc(H*W*sizeof(unsigned char));
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
		for(y = 1; y < H-1; y++)
			for(x = 1; x < W-1; x++)
				pnAccOut[y*W + x] = (unsigned char)((double)pnAccBl[y*W + x] / vmax * 255);

		if (flags&BDL_PGRAD_SAVEACC)
		{
			cname = (char*)malloc(256 * sizeof(char));
			strcpy(cname, name);
			if (CP == NULL)
				strcat(cname, "_acc1.bmp");
			else
				strcat(cname, "_acc2.bmp");
			CreateBmp8(cname, W, H, pnAccOut, 4);
			free(cname);
		}
		free(pnAccOut);
	}
	double vrmax = 0;
	vmax = rmax = 0;
	if(flags&BDL_PGRAD_CALCRADIUS)
	{
		if (CP != NULL)
		{
			rmn = (int)sqrt((double)rmn);
			rmx = (int)sqrt((double)rmx);
		}
		else
		{
			rmn = H + W;
			rmx = 0;
		}
		double* radHist = (double*)malloc((H+W)*sizeof(double));
		for(x = 0; x < H+W; x++)
			radHist[x] = 0.0;
		for(cGp=0;cGp<nGp;cGp++)
		{
			dx1 = (double)(pnEdgePnts[cGp].x - xmax);
			dy1 = (double)(pnEdgePnts[cGp].y - ymax);
			r1 = sqrt( dx1*dx1 + dy1*dy1 );
			if(r1 >= MINRAD & (r1 < rmn | r1 > rmx))
				radHist[(int)(r1+.5)]++;
		}
//		FILE* rhout = fopen("output.txt", "wt");
		for(x = 3; x < 398; x++)
			radHist[x] = 0.2 * radHist[x-2] + 0.2 * radHist[x-1] + 0.4 * radHist[x] + 0.2 * radHist[x+1] + 0.2 * radHist[x+2];
		for(x = 3; x < 398; x++){
			radHist[x] = 0.2 * radHist[x-2] + 0.2 * radHist[x-1] + 0.4 * radHist[x] + 0.2 * radHist[x+1] + 0.2 * radHist[x+2];
//			fprintf(rhout, "%1.2f  ", radHist[x]);
		}
//		fprintf(rhout, "\n");
//		fclose(rhout);

		if (CP != NULL)
		{
			for (x = 1; x < 400; x++)
				radHist[x] = radHist[x]/x;
		}
		for (x=1; x<400; x++)
			if (vrmax < radHist[x])
				vrmax = radHist[rmax = x];	
		free(radHist);
	}
	nBp = 0;
	vmax = 0;
/*	if (CP != NULL)
	{
		for (y = 1; y < H-1; y++)
			for (x = 1; x < W-1; x++)
			{
				cPt = y * W + x;
				gval = pgx_calc[cPt] * pgx_calc[cPt] + pgy_calc[cPt] * pgy_calc[cPt];
				val = (x-xmax) * (x-xmax) + (y-ymax) * (y-ymax);
				scalar = (x-xmax) * pgx_calc[cPt] + (y - ymax) * pgy_calc[cPt];

				if (gval > Tg && (double)(scalar*scalar) > 0.8 * (double)(val * gval)  && val <= (rmax+10)*(rmax+10) && val >= (rmax-10)*(rmax-10))
				{
					if (gval > vmax)
						vmax = gval;
					pnBoundary[nBp].x = x;
					pnBoundary[nBp].y = y;
					pnBoundary[nBp].val = 1;
					nBp++;
				}
			}
		printf("%d \n", vmax);
	}
	*/
	free(pnAccBl);
	free(pvBuf);

	CI->xc = xmax;
	CI->yc = ymax;
	CI->r = rmax;
	return 0;
}