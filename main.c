#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <windows.h>
#include <conio.h>
#include <time.h>
#include "header.h"
#include "errorcodes.h"
#include "ipl.h"
#include "unitypes.h"
#include "imio.h"
#include "bdl.h"

//#include "bpl.h"
//#include "bdl.h"
//#include "dbg.h"

unsigned char palitra[256][4];

int main(int argc, char** argv) 
{
	clock_t time;
	unsigned char *imgInput, *imgOutput;
	char *imgName, *imgOut, *name, *cname;
	char tmp[MAX_FILENAME];
	int i,k,nImg,res,H,W, ec,erc, e, erp, eri, angle, nref, flags;
	int errChist[100], errCRhist[100], errhist[100], errRPhist[100], errRIhist[100];
	CInfo pup, iri, CI, CP, rCP, cen;
	CErr err;
	
	FILE *in = fopen("../params_res.txt", "r");
	FILE *log = fopen("../data/res/0629/log.txt", "wb");
	FILE *plot = fopen("../data/res/0629/errors.txt", "wb");
	FILE *imgFile;
	double sum = 0.0, abs_err=.0, abs_errC = .0, rel_err = .0, rel_errC = .0;
	
	nImg = 1;//2000;
	memset(errhist,0,100*sizeof(int));
	memset(errChist,0,100*sizeof(int));
	memset(errCRhist,0,100*sizeof(int));
	memset(errRPhist,0,100*sizeof(int));
	memset(errRIhist,0,100*sizeof(int));
	name = (char*)malloc(256*sizeof(char));

	if (argc < 2)
	{
		printf("Warning: Markup file not stated. Using standart path '../params_res.txt'\n");
		if(!freopen("../params_res.txt","r", stdin))
		{
			printf("Error: Markup file doesn't exist.\n");
			return 1;
		}
	}
	else
	{
		if(argc<3)
			printf("Warning: Number of images to process is not stated. Processing the default number is %d.\n", nImg);
		else
			nImg = atoi(argv[argc-1]);
		freopen(argv[argc-2], "r", stdin);
		printf("Processing %d images from the markup file %s...\n", nImg,argv[argc-2]);
	}
	//340 338 207 549
	//Errors:
	//97 135 287 529 534 570 620 929 1716 489
	for(i = 0; i < 2167; i++)
		fgets(tmp,MAX_FILENAME,in);
	
	nref = 0;
	for(k = 0; k < nImg; k++)
	{
		printf("%d:\n", k);

		imgName = (char*)malloc(MAX_FILENAME);
		strcpy(imgName,"../data/");
		imgOut = (char*)malloc(MAX_FILENAME);
		strcpy(imgOut,"../data/res/0629/");
		
		//Read the expert markup data
		if((res = getExpertMarkup(in,&pup,&iri,imgName,imgOut))!=0)
		{
			printf("Error: Expert markup extraction failed\n");
			continue;
		}
		//else
			//fprintf(log, "Expert:\n(%d ; %d) r = %d\n(%d ; %d) r = %d\n", pup.xc, pup.yc, pup.r, iri.xc, iri.yc, iri.r);
		
		//Open an image file
		printf("%s\n", imgName);
		//fprintf(log, "%s\n", imgName);
		imgFile = fopen(imgName,"rb");
		if (imgFile == NULL){
			printf("Error: Cannot open the image file.\n");
			k--;
			nImg--;
			continue;
		}

		//Read the image data
		if ((res = readBmp8(imgFile, &imgInput,&H,&W))!=0)
		{
			printf("Error: Can't read the image data.\n");
			k--;
			nImg--;
			continue;
		}
		imgOutput = (unsigned char*)malloc(H*W*sizeof(unsigned char));
		memcpy(imgOutput, imgInput, H*W*sizeof(unsigned char));
		//Reduce the file extension
		memcpy(name,imgOut,strlen(imgOut)-4);
		name[strlen(imgOut)-4] = 0;
		fprintf(log, "name: %s\n", name);
		time = clock();
		
		//Iris segmentation function
		angle = 120;
		flags = BDL_PGRAD_CALCRADIUS|BDL_PUPREF_SAVECROP|BDL_PUPREF_SAVEPOLAR|BDL_PUPREF_SAVECROP;
		if ((res = ibd_IrisSegmentation(&CI,&CP,&rCP, name, imgInput,imgOutput,H,W, angle, log, flags)) != 0)
			printf("Error: Iris segmentation fail.\n");
		else
		{
			time = clock() - time;
		}
		//Axis convertation ???
		CI.yc = H - CI.yc;
		CP.yc = H - CP.yc;
		rCP.yc = H - rCP.yc;


		//Print the error output
		if ((res = calcError(&err, &CP, &CI, &pup, &iri)) == 0)
		{
			fprintf(log, "%f %f %f %f\n", err.sum, err.ctr, err.pup, err.iri);
			printf("No refinement: %f %f %f %f\n", err.sum, err.ctr, err.pup, err.iri);
			abs_err += err.sum * iri.r;
			abs_errC += err.ctr * iri.r;
		}
		e = (int)(err.sum * 100);
		ec = (int)(err.ctr * 100);
		erp = (int)(err.pup * 100);
		eri = (int)(err.iri * 100);
		if ((res = calcError(&err, &rCP, &CI, &pup, &iri)) == 0)
		{
			erc = (int)(err.ctr * 100);
			fprintf(log, "%f %f %f %f\n", err.sum, err.ctr, err.pup, err.iri);
			rel_err += err.sum;
			rel_errC += err.ctr;
			
			
			printf("Refinement: %f %f %f %f\n", err.sum, err.ctr, err.pup, err.iri);

		}
		CI.yc = H - CI.yc;
		CP.yc = H - CP.yc;
		rCP.yc = H - rCP.yc;

		if (erc >= ec)
			nref++;
		
		if (e < 100)
			errhist[e]++;
		printf("erc=%d\n", erc);
		if (ec >= 0 || erc > 20)
		{
			DRAW_CircleInGray(imgOutput, W, H, CI.xc, CI.yc, CI.r, 255);
			SaveBmp8(name, "_pupil_refined.bmp", W, H, imgOutput, 4);
		}


		if (ec < 100)
		{
			errChist[ec]++;
		}

		if (erc < 100)
			errCRhist[erc]++;
		
		if (erp < 100)
			errRPhist[erp]++;
		if (eri < 100)
			errRIhist[eri]++;
		
		printf("R_I = %d, expert: %d\n", CI.r, iri.r);

		printf("time: %f seconds.\n\n", (double)time/CLOCKS_PER_SEC);
		sum += (double)time/CLOCKS_PER_SEC;

		free(imgOut);
		free(imgOutput);
		free(imgInput);
		free(imgName);	
		fclose(imgFile);
	}

	fprintf(plot, "%d images\nSummed error:\n", nImg);
	for(i = 1; i < 100; i++)
	{
		errCRhist[i] += errCRhist[i-1];
		errChist[i] += errChist[i-1];
		errhist[i] += errhist[i-1];
		errRPhist[i] += errRPhist[i-1];
		errRIhist[i] += errRIhist[i-1];
		fprintf(plot, "%d ", errhist[i-1]);
	}
	fprintf(plot, "%d", errhist[99]);
	fprintf(plot, "\nOriginal center error:\n");
	for(i = 1; i < 100; i++)
		fprintf(plot, "%d ", errChist[i]);
	fprintf(plot, "\nCenter error after refinement:\n");
	for(i = 1; i < 100; i++)
		fprintf(plot, "%d ", errCRhist[i]);
	fprintf(plot, "\nPupil radius error(base):\n");
	for(i = 1; i < 100; i++)
		fprintf(plot, "%d ", errRPhist[i]);
	fprintf(plot, "\nIris radius error(base):\n");
	for(i = 1; i < 100; i++)
		fprintf(plot, "%d ", errRIhist[i]);
	fprintf(plot, "nrefined: %d\n", nref);
	fprintf(plot, "\nabs: %f, abs_c: %f, rel: %f, rel_c:%f\n", abs_err/nImg,abs_errC/nImg, rel_err/nImg, rel_errC/nImg);
	printf("%f %f\n%d image(s).\nTotal: %f\nAverage time: %f\n", 100 * (double)(errChist[9])/k, 100 * (double)(errChist[19])/k,k,sum, sum/k);
	printf("Summed error histogram: %f %f\n ", 100 * (double)(errhist[9])/k, 100 * (double)(errhist[19])/k);
	free(name);
	fclose(in);
	fclose(plot);
	fclose(log);
	//fclose(out);
	return 0;
}

