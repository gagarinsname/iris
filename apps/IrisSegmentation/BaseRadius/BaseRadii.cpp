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
#include <opencv2\opencv.hpp>

#include "GlobalParams.h"
#include "BaseRadii.h"




RESULT_CODE getMarkupData(SMarkup* cMarkupData, FILE* markupFile)
{
    int i, res;
    char* tmp = (char*)malloc(2 * MAX_FILENAME);
    if (!fgets(tmp, MAX_FILENAME, markupFile))
    {
        fprintf(stderr, "[ERROR] getMarkupData() cannot read another line from markup file.");
        return ERROR_WRONG_INPUT;
    }

    fscanf(markupFile, "%s", tmp);
    strcpy(cMarkupData->filename, tmp);

    fscanf(markupFile, "%s", tmp);
    cMarkupData->pupil.xc = atoi(tmp);
    cMarkupData->iris.xc = atoi(tmp);
    fscanf(markupFile, "%s", tmp);
    cMarkupData->pupil.yc = atoi(tmp);
    cMarkupData->iris.yc = atoi(tmp);

    if (cMarkupData->pupil.xc == 0 && cMarkupData->pupil.yc == 0)
    {
        fprintf(stderr, "[ERROR] getMarkupData() null line read.");
        return ERROR_WRONG_INPUT;
    }

    // skipping unnecessary input
    for (i = 0; i < 9; i++)
        fscanf(markupFile, "%s", tmp);

    fscanf(markupFile, "%s", tmp);
    cMarkupData->BazRad.pu_xc = atoi(tmp);
    fscanf(markupFile, "%s", tmp);
    cMarkupData->BazRad.pu_yc = atoi(tmp);
    fscanf(markupFile, "%s", tmp);
    cMarkupData->BazRad.pu_r = atoi(tmp);
    fscanf(markupFile, "%s", tmp);
    cMarkupData->BazRad.pu_q = atoi(tmp);

    fscanf(markupFile, "%s", tmp);
    cMarkupData->BazRad.ir_xc = atoi(tmp);
    fscanf(markupFile, "%s", tmp);
    cMarkupData->BazRad.ir_yc = atoi(tmp);
    fscanf(markupFile, "%s", tmp);
    cMarkupData->BazRad.ir_r = atoi(tmp);
    fscanf(markupFile, "%s", tmp);
    cMarkupData->BazRad.ir_q = atoi(tmp);
    free(tmp);
    return ERROR_OK;
}

void writeProjectionFeatureEntry(FILE* featureFile, SMarkup* cMarkup, int* pProjR, int* pProjL, int maxRadius)
{
    // writing file name, CenX, CenY and projection length
    fprintf(featureFile, "%s,%d,%d,%d,", cMarkup->filename, cMarkup->pupil.xc, cMarkup->pupil.yc, maxRadius);
    //writing features
    for (int i = 0; i < maxRadius; ++i)
        fprintf(featureFile, "%d,", pProjR[i]);
    for (int i = 0; i < maxRadius; ++i)
        fprintf(featureFile, "%d,", pProjL[i]);
    // writing targets
    fprintf(featureFile, "%d,%d\n", cMarkup->BazRad.pu_r, cMarkup->BazRad.ir_r);
}


RESULT_CODE GenerateProjectionFeatures(char* srcFolder, char* srcMarkupFile, char* dstFile)
{
    if (NULL == srcMarkupFile || NULL == dstFile)
    {
        fprintf(stderr, "[ERROR] GenerateProjectionFeatures() wrong input parameters.\n");
        return ERROR_NULL_POINTER;
    }

    FILE *markupFile = NULL, *featureFile = NULL;
    RESULT_CODE res = 0;
    char tmp[MAX_FILENAME];
    
    if ((markupFile = fopen(srcMarkupFile, "r")) == NULL)
    {
        fprintf(stderr, "[ERROR] GenerateProjectionFeatures() wrong filename for markup file.\n");
        return ERROR_WRONG_INPUT;
    }
    
    if ((featureFile = fopen(dstFile, "w")) == NULL)
    {
        fprintf(stderr, "[ERROR] GenerateProjectionFeatures() wrong filename for destination feature file.\n");
        return ERROR_WRONG_INPUT;
    }

    fgets(tmp, MAX_FILENAME, markupFile);
    char *imagePath = (char*)malloc(MAX_FILENAME * sizeof(char));
    
    SMarkup cMarkupData;
    int i = 0;
    res = 0;
    while((res = getMarkupData(&cMarkupData, markupFile)) == ERROR_OK)
    {
        memset(imagePath, 0, MAX_FILENAME * sizeof(char));
        strcpy(imagePath, srcFolder);

        // fprintf(stderr, "[ERROR]: getMarkupData markup extraction fail.\n");
        // return ERROR_WRONG_INPUT;
        
        strcat(imagePath, cMarkupData.filename);
        printf("image path %s\n", imagePath);
        //Read the image data
        cv::Mat srcImage = cv::imread(imagePath, CV_8U);
        int H = srcImage.rows, W = srcImage.cols;
        printf("Markup entry %d: %s\n", i,cMarkupData.filename);
        printf("CenX: %d CenY: %d\n", cMarkupData.pupil.xc, cMarkupData.pupil.yc);

        int *pProjR = (int*)malloc(W * sizeof(int));
        int *pProjL = (int*)malloc(W * sizeof(int));
        int outRadius = W / 2;
        int buflen = 0;
        if ((res = IPL_PROJ_FindHoughDonatorProjection7(
            pProjR, pProjL,
            (const uint8*)srcImage.data, W, H, cMarkupData.pupil.xc, cMarkupData.pupil.yc,
            20, &outRadius, 0,// dense analysis
            NULL, &buflen)))
        {
            fprintf(stderr, "[ERROR]: buffer size calculation failed.\n");
            return res;
        }
        printf("ProjectionCalc(): buf len: %d\n", buflen);
        void* buf = malloc(buflen);
        if ((res = IPL_PROJ_FindHoughDonatorProjection7(
            pProjR, pProjL,
            (const uint8*)srcImage.data, W, H, cMarkupData.pupil.xc, cMarkupData.pupil.yc,
            20, &outRadius, 0,// dense analysis
            buf, &buflen)))
        {
            fprintf(stderr, "[ERROR]: projection computation fails.\n");
            return res;
        }
        printf("ProjectionCalc(): out radius: %d", outRadius);
        writeProjectionFeatureEntry(featureFile, &cMarkupData, pProjR, pProjL, outRadius);    

            // cv::namedWindow("MarkupImage", CV_WINDOW_AUTOSIZE); //create a window with the name "MyWindow"
            // cv::imshow("MyWindow", srcImage); //display the image which is stored in the 'img' in the "MyWindow" window

            // cv::waitKey(0); //wait infinite time for a keypress

            // cv::destroyWindow("MarkupImage"); //destroy the window with the name, "MyWindow"
        printf("\n");
        i++;
    }
    printf("processing result: %d\n", res);

    return ERROR_OK;
}