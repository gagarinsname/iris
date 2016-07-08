#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <windows.h>
#include <conio.h>
#include "header.h"
#include "imio.h"
#include "stddefs.h"
#define roundf(x) floor(x + 0.5f)


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

	*img = (unsigned char*)malloc(*H * *W * sizeof(unsigned char));
	
	for(i = 0; i < *H; i++)
		for (j = 0; j < *W; j++)
		{
			unsigned char tmp;
			fread(&tmp, 1, 1, imgFile);
			(*img)[i * *W + j]=(unsigned int)tmp;
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
	for (i = 1; i < 256; i++)
	{
		Palette[i].rgbBlue = Palette[i-1].rgbBlue + 1;
		Palette[i].rgbGreen = Palette[i-1].rgbGreen + 1;
		Palette[i].rgbRed = Palette[i-1].rgbRed + 1;
	}
	WriteFile (hFile, Palette, 256 * sizeof (RGBQUAD), &RW, NULL);
	for (i = 0; i < Height; i++)
	{
		for (j = 0; j < Width; j++)
		{
			WriteFile (hFile, &map[i * Width + j], sizeof(color), &RW, NULL);
		}
		// Выровняем по границе
		WriteFile (hFile, Palette, (3 * Width) % 4, &RW, NULL);
	}
	CloseHandle(hFile);
	return 0;
}

int SaveBmp8 (char *fname, char* label, int Width,int Height, unsigned char* map, BYTE color)
{
	int res = 0;
	char* cname;
	cname = (char*)malloc(256 * sizeof(char));
	strcpy(cname, fname);
	strcat(cname, label);
	if ((res = CreateBmp8(cname, Width, Height, map, color)) != 0)
		printf("Error: Can not create .bmp file %s.\n", cname);
	free(cname);
}

 
void DRAW_Line(int* mat,int H, int W, int x1, int y1, int x2, int y2)
{
	//Костыль
	int dx,dy,lengthX,lengthY,length;
    dx = (x2 - x1 >= 0 ? 1 : -1);
    dy = (y2 - y1 >= 0 ? 1 : -1); 
    lengthX = abs(x2 - x1);
    lengthY = abs(y2 - y1);
	length = max(lengthX, lengthY);
 
	  if (x1 < 0 || x2 < 0 || y1 < 0 || y2 < 0 || x1 > W || x2 > W || y1 > H | y2 > H)
		return;


      if (length == 0)
      {
            mat[y1 * W + x1]++;
      }
 
      if (lengthY <= lengthX)
      {
            // Начальные значения
            int x = x1;
            int y = y1;
            int d = -lengthX;
 
            // Основной цикл
            length++;
            while(length--)
            {
				  mat[y * W + x]++;
                  x += dx;
                  d += 2 * lengthY;
                  if (d > 0) {
                        d -= 2 * lengthX;
                        y += dy;
                  }
            }
      }
      else
      {
            // Начальные значения
            int x = x1;
            int y = y1;
            int d = - lengthY;
 
            // Основной цикл
            length++;
            while(length--)
            {
                  mat[y* W + x]++;
                  y += dy;
                  d += 2 * lengthX;
                  if (d > 0) {
                        d -= 2 * lengthY;
                        x += dx;
                  }
            }
      }
}

// draw the circle of given color in grayscale image
void DRAW_CircleInGray(unsigned char* im, int W, int H, int xc, int yc, int r, unsigned char color)
{
	double step,phi,t;
	int x,y;
	for (x = xc-5; x < xc+5; x++)
	{
		if ((x>=0)&&(yc>=0)&&(x<W)&&(yc<H))
			im[yc*W+x] = color;
		if ((x>=0)&&(yc+1>=0)&&(x<W)&&(yc+1<H))
			im[(yc+1)*W+x] = color;
		if ((x>=0)&&(yc-1>=0)&&(x<W)&&(yc-1<H))
			im[(yc-1)*W+x] = color;
	}
	for (y = yc-5; y < yc+5; y++)
	{
		if ((xc>=0)&&(y>=0)&&(xc<W)&&(y<H))
			im[y*W+xc] = color;
		if ((xc+1>=0)&&(y>=0)&&(xc+1<W)&&(y<H))
			im[y*W+xc+1] = color;
		if ((xc-1>=0)&&(y>=0)&&(xc-1<W)&&(y<H))
			im[y*W+xc-1] = color;
	}
    step = PI/(4*r);
    for (phi=0.;phi<2*PI;phi+=step)
    {
      t = cos(phi)*r;
      if (t>0.)
        x = xc+(int)(t+.5);
      else
        x = xc+(int)(t-.5);
      t = sin(phi)*r;
      if (t>0.)
        y = yc+(int)(t+.5);
      else
        y = yc+(int)(t-.5);
      if ((x>=0)&&(y>=0)&&(x<W)&&(y<H))
        im[y*W+x] = color;
    }
}


void DRAW_SequenceInGray(unsigned char* im, int W, int H, int* xP, int* yP, int N, unsigned char color)
{
	int i;
    for (i=0;i < N; i++)
    {
      if ((xP[i]>=0)&&(yP[i]>=0)&&(xP[i]<W)&&(yP[i]<H))
        im[yP[i]*W+xP[i]] = color;
    }
}

void Dilate3x3Cross(unsigned char* dst, const unsigned char* src, int W, int H)
{
  int i,j,m,q;

  // upper boundary
  memset(dst,0,W);
  // body
  for (i=1;i<H-1;i++)
  {
    // left
    dst[i*W] = 0;
    // body
    for (j=1;j<W-1;j++)
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
