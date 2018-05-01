/*-----------------------------------------------------------------------*/
/*                                                                       */
/*      Created:  11 November 2017                                       */
/*      Purpose:  Iris library API                                       */
/*      Authors:  Iurii Efimov                                           */
/*                                                                       */
/*-----------------------------------------------------------------------*/

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "GlobalParams.h"
#include "BaseRadii.h"


RESULT_CODE myDetectBaseRadii(
    SPupilInfo* psPI,         // OUT: pupil data
    SIrisInfo* psII,          // OUT: iris data
    SCenterInfo* psCI,  // IN:  center data
    unsigned char* im,  // IN:  source image
    int W,                   // IN:  image size
    int H,                   //
    int* kernel3x3,           // edge detection kernel: kernel3x3[] = { -3, 0, 3, -10, 0, 10, -3, 0, 3 };
    int angMin,               // minimal projection estimation angle
    int angMax,               // maximal projection estimation angle
    int minR,                 // minimal radius
    int maxR,                 // maximal radius
    int thrHoriz,             // threshold for horizontal derivative
    float thrDot,             // threshold for circular shape: 0.8 - 1.0
    int* pProjR,              // circular projection - right side nums
    int* pProjL               // circular projection - left side nums
)
{
    int res = ERROR_OK;
    // int maxR = MIA_min(MIA_max(MIA_min(psCI->xc, W - psCI->xc), MIA_min(psCI->yc, H - psCI->yc)), H / 2);

    if ((res = FindHoughProjection(pProjR, pProjL, angMin, angMax, kernel3x3, thrHoriz,thrDot,minR,maxR, im, W, H, psCI->xc, psCI->yc)))
    {
        fprintf(stderr, "[ERROR]: Projection feature generation failed.\n");
        return res;
    }
    memset(pProjL, 0, MINIMAL_PUPIL_RADIUS * sizeof(int));
    memset(pProjR, 0, MINIMAL_PUPIL_RADIUS * sizeof(int));

    int mode = IVIR_ARRPOXIR_PRPRPR | IVIR_ARRPOXIR_USEPROJ; //  IVIR_ARRPOXIR_LEFEYE | IVIR_ARRPOXIR_PRPRPR;
    if ((res = IVIR_PrPu(
        psPI, psII, psCI, im, W, H,
        mode, -1, -1, pProjR, pProjL, NULL)))
    {
        fprintf(stderr, "[ERROR]: base radii estimations failed.\n");
    }
    return res;
}


int DetectBaseRadii(
    SPupilInfo* psPI,         // OUT: pupil data
    SIrisInfo* psII,          // OUT: iris data
    const SCenterInfo* psCI,  // IN:  center data
    const unsigned char* im,  // IN:  source image
    int xs,                   // IN:  image size
    int ys,                   //
    int* pProjR,               // circular projection - right side nums
    int* pProjL               // circular projection - left side nums
                              // int mode // mask: 1- draw projs to buf, 2- use projs from buf
)
{
    RESULT_CODE res;
    int outRadius = xs / 2;
    int buflen = 0;

    if ((res = IPL_PROJ_FindHoughDonatorProjection7(
        pProjR, pProjL,
        im, xs, ys, psCI->xc, psCI->yc,
        MINIMAL_PUPIL_RADIUS, &outRadius, 0,// dense analysis
        NULL, &buflen)))
    {
        fprintf(stderr, "[ERROR]: buffer size calculation failed.\n");
        return res;
    }
    void* buf = malloc(buflen);
    if ((res = IPL_PROJ_FindHoughDonatorProjection7(
        pProjR, pProjL,
        im, xs, ys, psCI->xc, psCI->yc,
        MINIMAL_PUPIL_RADIUS, &outRadius, 0,// dense analysis
        buf, &buflen)))
    {
        free(buf);
        fprintf(stderr, "[ERROR]: projection computation fails.\n");
        return res;
    }

    int mode = IVIR_ARRPOXIR_PRPRPR; //  IVIR_ARRPOXIR_LEFEYE | IVIR_ARRPOXIR_PRPRPR;
    if ((res = IVIR_PrPu(
        psPI, psII, psCI, im, xs, ys,
        mode, -1, -1, pProjR, pProjL, NULL)))
    {
        free(buf);
        fprintf(stderr, "[ERROR]: base radii estimations failed.\n");
        return res;
    }
    free(buf);
    return ERROR_OK;
}
