#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <windows.h>
#include <conio.h>
#include "header.h"
#include "imio.h"
#include "stddefs.h"


int readBmp8(FILE* imgFile, unsigned char** img,int* H, int* W)
{
	BITMAPFILEHEADER bfh;
    BITMAPINFOHEADER bih;
	//*img = NULL;
	int i,j;
	if(!(fread(&bfh, sizeof(bfh), 1, imgFile))){
			printf("Error: Cannot read the header\n");
			return(-1);
	}

	if(!(fread(&bih, sizeof(bih), 1, imgFile))){
		printf("Error: Cannot read header\n");
		return(-1);
	}
	//printf("%d x %d, %d bit(s),\n", bih.biWidth, bih.biHeight, bih.biBitCount);

	if (bih.biBitCount != 8)
	{
		printf("Error: Wrong image bit number.\n");
		return -1;
	}
	*H = bih.biHeight;
	*W = bih.biWidth;

	fseek(imgFile, bfh.bfOffBits, SEEK_SET);

	*img = (uint8*)malloc(*H * *W * sizeof(uint8));
	
	for(i = *H-1; i > -1; --i)
		for (j = 0; j < *W; ++j)
		{
			unsigned char tmp;
			fread(&tmp, 1, 1, imgFile);
			(*img)[i * *W + j]=(uint8)tmp;
		}

	return 0;
}

int getExpertMarkup(FILE* fin, CInfo* pup,CInfo* iri,char* name,char* out)
{
	int i;
	char* tmp = (char*)malloc(10*MAX_FILENAME);
	fgets(tmp,MAX_FILENAME,fin);
	fscanf(fin, "%s", tmp);
	strcat(name,tmp);
	strcat(out,tmp);
	//printf("%s\n", name);
	for(i = 0; i < 11; i++)
	{
		fscanf(fin, "%s", tmp);
	}
	fscanf(fin, "%s", tmp);
	pup->xc = atoi(tmp);
	fscanf(fin, "%s", tmp);
	pup->yc = atoi(tmp);
	fscanf(fin, "%s", tmp);
	pup->r = atoi(tmp);
	fscanf(fin, "%s", tmp);
	fscanf(fin, "%s", tmp);
	iri->xc = atoi(tmp);
	fscanf(fin, "%s", tmp);
	iri->yc = atoi(tmp);
	fscanf(fin, "%s", tmp);
	iri->r = atoi(tmp);
//	printf("Expert:\n(%d ; %d) r = %d\n(%d ; %d) r = %d\n", pup->xc, pup->yc, pup->r, iri->xc, iri->yc, iri->r);
	free(tmp);
	return 0;
}

int calcError(CErr* err, CInfo* myPup, CInfo* myIri, CInfo* pup, CInfo* iri)
{
	double d1, d2;
	if (myIri->r < myPup->r)
	{
		int16 tmp;
		tmp = myIri->r;
		myIri->r = myPup->r;
		myPup->r = tmp;
		tmp = myIri->xc;
		myIri->xc = myPup->xc;
		myPup->xc = tmp;
		tmp = myIri->yc;
		myIri->yc = myPup->yc;
		myPup->yc = tmp;
	}
	printf("pup: (%d , %d , %d)\n iri: (%d , %d , %d)\n", myPup->xc, myPup->yc, myPup->r, myIri->xc, myIri->yc, myIri->r);
	err->iri = (double)(abs(myIri->r - iri->r));
	err->iri /= iri->r;
	err->pup = (double)(abs(myPup->r - pup->r));
	err->pup /= iri->r;
	
	d1 = sqrt((double)((myPup->xc - pup->xc)*(myPup->xc - pup->xc) + (myPup->yc - pup->yc)*(myPup->yc - pup->yc)));// + (myIri->xc - iri->xc)*(myIri->xc - iri->xc) + (myIri->yc - iri->yc)*(myIri->yc - iri->yc));
	//d2 = sqrt((double)();
	err->ctr = d1;
	err->ctr /= iri->r;
	err->sum = err->ctr + err->iri + err->pup;

	return 0;
}

int CreateBmp8 (char *fname, int Width,int Height, unsigned char* map, BYTE color)
{
	HANDLE hFile;
	DWORD RW;
	int i, j;

	// Объявим нужные структуры
	BITMAPFILEHEADER bfh;
	BITMAPINFOHEADER bih;
	RGBQUAD Palette [256];// Палитра

	// Заполним их
	memset (&bfh, 0, sizeof(bfh));
	bfh.bfType = 0x4D42;									// Обозначим, что это bmp 'BM'
	bfh.bfOffBits = sizeof(bfh) + sizeof(bih) + 1024;		// Палитра занимает 1Kb, но мы его испоьзовать не будем
	bfh.bfSize = bfh.bfOffBits + 
				 sizeof(color) * Width * Height + 
				 Height * ((3*Width) % 4);					// Посчитаем размер конечного файла

	memset (&bih, 0, sizeof(bih));
	bih.biSize = sizeof(bih);								// Так положено
	bih.biBitCount = 8;										// 16 бит на пиксель
	bih.biCompression = BI_RGB;								// Без сжатия
	bih.biHeight = Height;
	bih.biWidth = Width;
	bih.biPlanes = 1;										// Должно быть 1
															// А остальные поля остаются 0
	hFile = CreateFile(fname, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, 0, NULL);
	if (hFile == INVALID_HANDLE_VALUE)
		return -1;

	// Запишем заголовки
	WriteFile (hFile, &bfh, sizeof (bfh), &RW, NULL);
	WriteFile (hFile, &bih, sizeof (bih), &RW, NULL);

	// Создадим и запишем палитру
	memset (&Palette[0], 0, sizeof (RGBQUAD));
	for (i = 1; i < 256; ++i)
	{
		Palette[i].rgbBlue = Palette[i-1].rgbBlue + 1;
		Palette[i].rgbGreen = Palette[i-1].rgbGreen + 1;
		Palette[i].rgbRed = Palette[i-1].rgbRed + 1;
	}
	WriteFile (hFile, Palette, 256 * sizeof (RGBQUAD), &RW, NULL);
	for (i = Height-1; i > -1; --i)
	{
		for (j = 0; j < Width; ++j)
		{
			WriteFile (hFile, &map[i * Width + j], sizeof(color), &RW, NULL);
		}
		// Выровняем по границе
		WriteFile (hFile, Palette, (3 * Width) % 4, &RW, NULL);
	}
	CloseHandle(hFile);
	return 0;
}

int SaveBmp8 (char *fname, char* label, int Width, int Height, unsigned char* map, BYTE color)
{
	int res = 0;
	char* cname;
	char* ext;
	cname = (char*)malloc(MAX_FILENAME * sizeof(char));
	ext = (char*)malloc(MAX_FILENAME * sizeof(char));
	strcpy(ext, fname + strlen(fname) - 4);
	memcpy(cname, fname, strlen(fname) - 4);
	cname[strlen(fname) - 4] = 0;
	 
	strcat(cname, label);
	strcat(cname, ext);

	if ((res = CreateBmp8(cname, Width, Height, map, color)) != 0)
	{
		fprintf(stderr, "[ERROR]: SaveBmp8() can not create .bmp file %s.\n", cname);
		return -1;
	}

	free(cname);
	free(ext);
	return res;
}



void Dilate3x3Cross(unsigned char* dst, const unsigned char* src, int W, int H)
{
  int i,j,m,q;

  // upper boundary
  memset(dst,0,W);
  // body
  for (i=1;i<H-1;++i)
  {
    // left
    dst[i*W] = 0;
    // body
    for (j=1;j<W-1;++j)
    {
      m =  src[(i-1)*W+(j  )];

      if (m<(q = src[(i  )*W+(j-1)]))
        m = q;
      if (m<(q = src[(i  )*W+(j  )]))
        m = q;
      if (m<(q = src[(i  )*W+(j+1)]))
        m = q;
      if (m<(q = src[(i+1)*W+(j  )]))
        m = q;
      dst[i*W+j] = (unsigned char)m;
    }
    // right
    dst[i*W+W-1] = 0;
  }
  // lower boundary
  memset(dst+(H-1)*W,0,W);
}
