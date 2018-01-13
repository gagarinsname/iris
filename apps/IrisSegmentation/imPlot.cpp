#include <math.h>
#include <stdio.h>
#include "imPlot.h"
#include "unitypes.h"

void DRAW_Line(int* mat, int H, int W, int x1, int y1, int x2, int y2)
{
	//Костыль
	int dx, dy, lengthX, lengthY, length;
	int longLength, shortLength;
	int x, y, d;

	dx = (x2 - x1 >= 0 ? 1 : -1);
	dy = (y2 - y1 >= 0 ? 1 : -1);
	lengthX = abs(x2 - x1);
	lengthY = abs(y2 - y1);
		
	length = fmax(lengthX, lengthY);

	if (length == 0)
		mat[y1 * W + x1]++;

	if (lengthY <= lengthX)
	{
		// Начальные значения
		x = x1;
		y = y1;
		d = -lengthX;

		// Основной цикл
		++length;
		while (length--)
		{
			if (x > -1 && x < W && y > -1 && y < H)
				++(mat[y * W + x]);

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
		x = x1;
		y = y1;
		d = -lengthY;

		// Основной цикл
		++length;
		while (--length)
		{
			if (x > -1 && x < W && y > -1 && y < H)
				++(mat[y* W + x]);

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
void DRAW_CircleInGray(
	uint8* im,
	int W,
	int H,
	SCircleData* sCircle,
	uint8 color)
{
	double step, phi, t;
	int x, y;
	for (x = sCircle->xc - 5; x < sCircle->xc + 5; ++x)
	{
		if ((x >= 0) && (sCircle->yc >= 0) && (x < W) && (sCircle->yc < H))
			im[sCircle->yc * W + x] = color;
		if ((x >= 0) && (sCircle->yc + 1 >= 0) && (x < W) && (sCircle->yc + 1 < H))
			im[(sCircle->yc + 1) * W + x] = color;
		if ((x >= 0) && (sCircle->yc - 1 >= 0) && (x < W) && (sCircle->yc - 1 < H))
			im[(sCircle->yc - 1) * W + x] = color;
	}
	for (y = sCircle->yc - 5; y < sCircle->yc + 5; ++y)
	{
		if ((sCircle->xc >= 0) && (y >= 0) && (sCircle->xc<W) && (y<H))
			im[y * W + sCircle->xc] = color;
		if ((sCircle->xc + 1 >= 0) && (y >= 0) && (sCircle->xc + 1 < W) && (y < H))
			im[y * W + sCircle->xc + 1] = color;
		if ((sCircle->xc - 1 >= 0) && (y >= 0) && (sCircle->xc - 1 < W) && (y < H))
			im[y * W + sCircle->xc - 1] = color;
	}
	step = PI / (4 * sCircle->r);
	for (phi = 0.; phi < 2 * PI; phi += step)
	{
		t = cos(phi)*sCircle->r;
		if (t>0.)
			x = sCircle->xc + (int)(t + .5);
		else
			x = sCircle->xc + (int)(t - .5);
		t = sin(phi) * sCircle->r;
		if (t > 0.)
			y = sCircle->yc + (int)(t + .5);
		else
			y = sCircle->yc + (int)(t - .5);
		
		if ((x >= 0) && (y >= 0) && (x < W) && (y < H))
			im[y * W + x] = color;
	}
}

void DRAW_2DSequenceInGray(uint8* im, int W, int H, sPoint* sSeq, int N, uint8 color)
{
	for (int i = 0; i < N; ++i)
	{
		if ((sSeq[i].x >= 0) && (sSeq[i].y >= 0) && (sSeq[i].x < W) && (sSeq[i].y < H))
			im[sSeq[i].y * W + sSeq[i].x] = color;
	}
}
