#include <stdio.h>
#include <memory.h>
#include <malloc.h>

#include "ImageProcessing.h"

// hist processing
RESULT_CODE IPL_HIST_CalcMedianInLine(
    unsigned char* line,
    const unsigned char* im,
    int xs,
    int ys,
    int h,
    int kew)
{
    int i, hist[256], vmed, dx, dy, dxb, nlef, nrig, v;

    // check arguments
    if (2 * kew + 1>ys)
        return ERROR_WRONG_INPUT;
    if (h<kew)
        h = kew;
    if (h + kew >= ys)
        h = ys - kew - 1;
    // initial set
    im += h*xs;
    memset(hist, 0, sizeof(hist));
    for (dx = -kew; dx<kew; dx++) // '<' rather than '<='
    {
        dxb = MIA_bound(dx, 0, xs - 1);
        for (dy = -kew; dy <= kew; dy++)
            hist[im[dy*xs + dxb]]++;
    }
    vmed = 0;
    nlef = 0;
    nrig = (2 * kew + 1)*(2 * kew) - hist[0];
    // main cycle
    for (i = 0; i<xs; i++)
    {
        // add points from right
        dxb = MIA_bound(i + kew, 0, xs - 1);
        for (dy = -kew; dy <= kew; dy++)
        {
            hist[v = im[dy*xs + dxb]]++;
            if (v>vmed)
                nrig++;
            else
                if (v<vmed)
                    nlef++;
        }
        // move median
        while (nrig>hist[vmed] + nlef)
        {
            nlef += hist[vmed];
            vmed++;
            nrig -= hist[vmed];
        }
        while (nlef>hist[vmed] + nrig)
        {
            nrig += hist[vmed];
            vmed--;
            nlef -= hist[vmed];
        }
        // save value
        line[i] = (unsigned char)vmed;
        // remove points from left
        dxb = MIA_bound(i - kew, 0, xs - 1);
        for (dy = -kew; dy <= kew; dy++)
        {
            hist[v = im[dy*xs + dxb]]--;
            if (v>vmed)
                nrig--;
            else
                if (v<vmed)
                    nlef--;
        }
    }
    return ERROR_OK;
}

void IPL_HIST_Blur(
    int* dst,         // OUT: blurred hist
    const int* val,   // IN:  source hist
    int len,          // IN:  length of the hist
    int hw)           // IN:  half-window
{
    int i, V;

    // beginning
    V = val[0] * (hw + 1);
    for (i = 1; i<hw; i++)
        V += val[i];
    for (i = 0; i<hw; i++)
    {
        V += val[i + hw];
        dst[i] = V;
        V -= val[0];
    }
    // body
    len -= hw;
    for (; i<len; i++)
    {
        V += val[i + hw];
        dst[i] = V;
        V -= val[i - hw];
    }
    // ending
    len += hw;
    for (; i<len; i++)
    {
        V += val[len - 1];
        dst[i] = V;
        V -= val[i - hw];
    }
}

void IPL_HIST_Blur_double(
    double* dst,         // OUT: blurred hist
    const double* val,   // IN:  source hist
    int len,          // IN:  length of the hist
    int hw)           // IN:  half-window
{
    int i;
    double V;

    // beginning
    V = val[0] * (hw + 1);
    for (i = 1; i<hw; i++)
        V += val[i];
    for (i = 0; i<hw; i++)
    {
        V += val[i + hw];
        dst[i] = V;
        V -= val[0];
    }
    // body
    len -= hw;
    for (; i<len; i++)
    {
        V += val[i + hw];
        dst[i] = V;
        V -= val[i - hw];
    }
    // ending
    len += hw;
    for (; i<len; i++)
    {
        V += val[len - 1];
        dst[i] = V;
        V -= val[i - hw];
    }
}


// make median of all array for byte
int IPL_HIST_mdn_pixE(
	const unsigned char* bufin, // input array
	int                  lstr)  // size bufin
{                             // Returns: median
	int i, half, S, anARM[256];

	memset(&anARM[0], 0, sizeof(anARM));
	half = (lstr / 2) | 1;
	for (i = 0; i<lstr; anARM[bufin[i++]]++);
	for (i = 0, S = 0; S<half; S += anARM[i++]);
	return i - 1;
}

/****************** IMAGE FILTERING *****************************************/
// Dilatation with 3x3 cross kernel
void Dilate3x3Cross(unsigned char* dst, const unsigned char* src, int W, int H)
{
	int i, j, m, q;

	// upper boundary
	memset(dst, 0, W);
	// body
	for (i = 1; i<H - 1; ++i)
	{
		// left
		dst[i*W] = 0;
		// body
		for (j = 1; j<W - 1; ++j)
		{
			m = src[(i - 1)*W + (j)];

			if (m<(q = src[(i)*W + (j - 1)]))
				m = q;
			if (m<(q = src[(i)*W + (j)]))
				m = q;
			if (m<(q = src[(i)*W + (j + 1)]))
				m = q;
			if (m<(q = src[(i + 1)*W + (j)]))
				m = q;
			dst[i*W + j] = (unsigned char)m;
		}
		// right
		dst[i*W + W - 1] = 0;
	}
	// lower boundary
	memset(dst + (H - 1)*W, 0, W);
}


#define BOOSTBLURFACTOR 90.0

// Blur an image with a gaussian filter.
RESULT_CODE IPL_FILT_GaussianSmooth_uint8(
    short *smoothedim,
    double *tempim,
    const unsigned char *image,
    int xs,
    int ys,
    double sigma)
{
    int r, c, rr, cc, center, i;
    double *kernel, dot, sum, fx;

    if ((smoothedim == NULL) || (tempim == NULL) || (image == NULL))
        return ERROR_NULL_POINTER;
    if ((xs <= 0) || (ys <= 0) || (sigma <= 0.))
        return ERROR_WRONG_INPUT;
    // Create a 1-dimensional Haussian smoothing kernel.
    sum = 0.0;
    center = (int)ceil(2.5 * sigma);
    if ((kernel = (double*)calloc(1 + 2 * center, sizeof(double))) == NULL)
        return ERROR_MEMORY;
    for (i = -center; i <= center; i++)
    {
        fx = (double)(pow(2.71828, -0.5*i*i / (sigma*sigma)) / (sigma*sqrt(6.2831853)));
        kernel[i + center] = fx;
        sum += fx;
    }
    for (i = -center; i <= center; i++)
        kernel[i + center] /= sum;
    // Blur in the x - direction.
    for (r = 0; r<ys; r++)
        for (c = 0; c<xs; c++)
        {
            dot = 0.0;
            sum = 0.0;
            for (cc = (-center); cc <= center; cc++)
                if (((c + cc) >= 0) && ((c + cc) < xs))
                {
                    dot += (double)image[r*xs + (c + cc)] * kernel[center + cc];
                    sum += kernel[center + cc];
                }
            tempim[r*xs + c] = dot / sum;
        }
    // Blur in the y - direction.
    for (c = 0; c<xs; c++)
        for (r = 0; r<ys; r++)
        {
            sum = 0.0;
            dot = 0.0;
            for (rr = (-center); rr <= center; rr++)
                if (((r + rr) >= 0) && ((r + rr) < ys))
                {
                    dot += tempim[(r + rr)*xs + c] * kernel[center + rr];
                    sum += kernel[center + rr];
                }
            smoothedim[r*xs + c] = (short)(dot*BOOSTBLURFACTOR / sum + 0.5);
        }
    free(kernel);
    return ERROR_OK;
}

// Blur an image with a gaussian filter
RESULT_CODE IPL_FILT_GaussianSmooth_ext_uint8(
    short *smoothedim,
    double *tempim,
    const unsigned char *image,
    int xs,
    int ys,
    double sigma_x,
    double sigma_y,
    double _BOOSTBLURFACTOR)
{
    int r, c, rr, cc, center, i;
    double *kernel, dot, sum, fx;

    if ((smoothedim == NULL) || (tempim == NULL) || (image == NULL))
        return ERROR_NULL_POINTER;
    if ((xs <= 0) || (ys <= 0) || (sigma_x <= 0.) || (sigma_y <= 0.))
        return ERROR_WRONG_INPUT;
    // Blur in the x - direction.
    // Create a 1-dimensional Haussian smoothing kernel.
    sum = 0.0;
    center = (int)ceil(2.5 * sigma_x);
    if ((kernel = (double*)calloc(1 + 2 * center, sizeof(double))) == NULL)
        return ERROR_MEMORY;
    for (i = -center; i <= center; i++)
    {
        fx = (double)(pow(2.71828, -0.5*i*i / (sigma_x*sigma_x)) / (sigma_x*sqrt(6.2831853)));
        kernel[i + center] = fx;
        sum += fx;
    }
    for (i = -center; i <= center; i++)
        kernel[i + center] /= sum;
    // Blur
    for (r = 0; r<ys; r++)
        for (c = 0; c<xs; c++)
        {
            dot = 0.0;
            sum = 0.0;
            for (cc = (-center); cc <= center; cc++)
                if (((c + cc) >= 0) && ((c + cc) < xs))
                {
                    dot += (double)image[r*xs + (c + cc)] * kernel[center + cc];
                    sum += kernel[center + cc];
                }
            tempim[r*xs + c] = dot / sum;
        }
    free(kernel);
    // Blur in the y - direction.
    // Create a 1-dimensional Haussian smoothing kernel.
    sum = 0.0;
    center = (int)ceil(2.5 * sigma_y);
    if ((kernel = (double*)calloc(1 + 2 * center, sizeof(double))) == NULL)
        return ERROR_MEMORY;
    for (i = -center; i <= center; i++)
    {
        fx = (double)(pow(2.71828, -0.5*i*i / (sigma_y*sigma_y)) / (sigma_y*sqrt(6.2831853)));
        kernel[i + center] = fx;
        sum += fx;
    }
    for (i = -center; i <= center; i++)
        kernel[i + center] /= sum;
    // Blur
    for (c = 0; c<xs; c++)
        for (r = 0; r<ys; r++)
        {
            sum = 0.0;
            dot = 0.0;
            for (rr = (-center); rr <= center; rr++)
                if (((r + rr) >= 0) && ((r + rr) < ys))
                {
                    dot += tempim[(r + rr)*xs + c] * kernel[center + rr];
                    sum += kernel[center + rr];
                }
            smoothedim[r*xs + c] = (short)(dot*_BOOSTBLURFACTOR / sum + 0.5);
        }
    free(kernel);
    return ERROR_OK;
}

// Blur an image with a gaussian filter.
RESULT_CODE IPL_FILT_GaussianSmooth_double(
	double *smoothedim,
	double *tempim,
	const double *image,
	int xs,
	int ys,
	double sigma)
{
	int r, c, rr, cc, center, i;
	double *kernel, dot, sum, fx;

	if ((smoothedim == NULL) || (tempim == NULL) || (image == NULL))
		return ERROR_NULL_POINTER;
	if ((xs <= 0) || (ys <= 0) || (sigma <= 0.))
		return ERROR_WRONG_INPUT;
	// Create a 1-dimensional gaussian smoothing kernel.
	sum = 0.0;
	center = (int)ceil(2.5 * sigma);
	if ((kernel = (double*)calloc(1 + 2 * center, sizeof(double))) == NULL)
		return ERROR_MEMORY;
	for (i = -center; i <= center; i++)
	{
		fx = (double)(pow(2.71828, -0.5*i*i / (sigma*sigma)) / (sigma*sqrt(6.2831853)));
		kernel[i + center] = fx;
		sum += fx;
	}
	for (i = -center; i <= center; i++)
		kernel[i + center] /= sum;
	// Blur in the x - direction.
	for (r = 0; r < ys; r++)
		for (c = 0; c < xs; c++)
		{
			dot = 0.0;
			sum = 0.0;
			for (cc = (-center); cc <= center; cc++)
				if (((c + cc) >= 0) && ((c + cc) < xs))
				{
					dot += image[r*xs + (c + cc)] * kernel[center + cc];
					sum += kernel[center + cc];
				}
			tempim[r*xs + c] = dot / sum;
		}
	// Blur in the y - direction.
	for (c = 0; c < xs; c++)
		for (r = 0; r < ys; r++)
		{
			sum = 0.0;
			dot = 0.0;
			for (rr = (-center); rr <= center; rr++)
				if (((r + rr) >= 0) && ((r + rr) < ys))
				{
					dot += tempim[(r + rr)*xs + c] * kernel[center + rr];
					sum += kernel[center + rr];
				}
			smoothedim[r*xs + c] = dot*BOOSTBLURFACTOR / sum;
		}
	free(kernel);
	return ERROR_OK;
}

// Sobel filtering
int IPL_Sobel(uint8* pucG, int16* pGx, int16* pGy, const uint8* img, int H, int W, int ksize, int *KernelKxK, int pad)
{
	if (NULL == pucG || NULL == pGx || NULL == pGy || NULL == img || NULL == KernelKxK)
		return ERROR_NULL_POINTER;
	if (H <= 0 || W <= 0 || ksize <= 3)
		return ERROR_WRONG_INPUT;
	if (pad < ksize)
		pad = ksize;

	return ERROR_OK;
}


/******************* GAUSSIAN FILTERING FROM MIA **********************/
void ownIPL_FillBorders1B(
	unsigned char* im,
	int xs,
	int ys)
{
	int i;

	// lef & rig
	for (i = 1; i<ys - 1; i++)
	{
		im[i*xs] = im[i*xs + 1];
		im[i*xs + xs - 1] = im[i*xs + xs - 2];
	}
	// top & bot
	memcpy(im, im + xs, xs);
	memcpy(im + (ys - 1)*xs, im + (ys - 2)*xs, xs);
}


void ownIPL_FillBorders1B_int(
	int* im,
	int xs,
	int ys)
{
	int i;

	// lef & rig
	for (i = 1; i<ys - 1; i++)
	{
		im[i*xs] = im[i*xs + 1];
		im[i*xs + xs - 1] = im[i*xs + xs - 2];
	}
	// top & bot
	memcpy(im, im + xs, xs * sizeof(im[0]));
	memcpy(im + (ys - 1)*xs, im + (ys - 2)*xs, xs * sizeof(im[0]));
}


void ownIPL_FillBordersXB(
	unsigned char* dst,
	int xs,
	int ys,
	int kehw,
	int kehh)
{
	int i, j;

	// sides
	for (i = kehh; i<ys - kehh; i++)
		for (j = 0; j<kehw; j++)
		{
			dst[i*xs + j] = dst[i*xs + kehw];
			dst[i*xs + (xs - kehw) + j] = dst[i*xs + (xs - kehw) - 1];
		}
	// top & bottom
	for (i = 0; i<kehh; i++)
	{
		memcpy(dst + i*xs, dst + kehh*xs, xs);
		memcpy(dst + (ys - kehh + i)*xs, dst + (ys - kehh - 1)*xs, xs);
	}
}

void IPL_FILT_HaussBlur3x3(
	unsigned char* dst,
	const unsigned char* src,
	int xs,
	int ys)
{
	int i, j;

	for (i = 1; i<ys - 1; i++)
	{
		for (j = 1; j<xs - 1; j++)
			dst[i*xs + j] = (unsigned char)(
			(1 * src[(i - 1)*xs + (j - 1)] + 2 * src[(i - 1)*xs + (j)] + 1 * src[(i - 1)*xs + (j + 1)] +
				2 * src[(i)*xs + (j - 1)] + 4 * src[(i)*xs + (j)] + 2 * src[(i)*xs + (j + 1)] +
				1 * src[(i + 1)*xs + (j - 1)] + 2 * src[(i + 1)*xs + (j)] + 1 * src[(i + 1)*xs + (j + 1)]) / 16);
	}
	ownIPL_FillBorders1B(dst, xs, ys);
}


void IPL_FILT_HaussBlur3x3_Int(
	int* dst,
	const int* src,
	int xs,
	int ys,
	int mode) // 0 - normalized, 1- no normalization
{
	int i;

	if (mode)
	{ // no normalization
	  // body
		dst += xs + 1;
		src += xs + 1;
		for (i = (ys - 2)*xs - 2; i>0; i--)
		{
			*dst = 1 * src[-xs - 1] + 2 * src[-xs] + 1 * src[-xs + 1] +
				2 * src[-1] + 4 * src[0] + 2 * src[1] +
				1 * src[xs - 1] + 2 * src[xs] + 1 * src[xs + 1];
			src++;
			dst++;
		}
		dst -= xs*ys - xs - 1;
	}
	else
	{ // normalized
	  // body
		dst += xs + 1;
		src += xs + 1;
		for (i = (ys - 2)*xs - 2; i>0; i--)
		{
			*dst = (1 * src[-xs - 1] + 2 * src[-xs] + 1 * src[-xs + 1] +
				2 * src[-1] + 4 * src[0] + 2 * src[1] +
				1 * src[xs - 1] + 2 * src[xs] + 1 * src[xs + 1]) / 16;
			src++;
			dst++;
		}
		dst -= xs*ys - xs - 1;
	}
	// borders
	ownIPL_FillBorders1B_int(dst, xs, ys);
}


void IPL_FILT_HaussBlur1x3(
	unsigned char* dst,
	const unsigned char* src,
	int xs,
	int ys)
{
	int i, j;

	for (i = 1; i<ys - 1; i++)
	{
		for (j = 0; j<xs - 0; j++)
			dst[i*xs + j] = (unsigned char)((52 * src[(i - 1)*xs + (j)] + 127 * src[(i)*xs + (j)] + 52 * src[(i + 1)*xs + (j)]) / 231);
		//dst[i*xs+j] = (int)( src[(i  )*xs+(j  )]);
	}
	ownIPL_FillBordersXB(dst, xs, ys, 0, 1);
}


void IPL_FILT_HaussBlur3x1(
	unsigned char* dst,
	const unsigned char* src,
	int xs,
	int ys)
{
	int i, j;

	for (i = 0; i<ys - 0; i++)
	{
		for (j = 1; j<xs - 1; j++)
			dst[i*xs + j] = (unsigned char)((52 * src[(i)*xs + (j - 1)] + 127 * src[(i)*xs + (j)] + 52 * src[(i)*xs + (j + 1)]) / 231);
		//dst[i*xs+j] = (int)( src[(i  )*xs+(j  )]);
	}
	ownIPL_FillBordersXB(dst, xs, ys, 1, 0);
}
