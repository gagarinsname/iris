#include <stdio.h>
#include <memory.h>
#include <malloc.h>

// #include "GlobalParams.h"
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


// Image filtering

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
    for (r = 0; r<ys; r++)
        for (c = 0; c<xs; c++)
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
            smoothedim[r*xs + c] = dot*BOOSTBLURFACTOR / sum;
        }
    free(kernel);
    return ERROR_OK;
}