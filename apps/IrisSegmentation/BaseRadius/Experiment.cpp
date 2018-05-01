#include <stdio.h>
#include <string.h>
#include <iostream>
#include "GlobalParams.h"
#include "BaseRadii.h"

#include <opencv2/opencv.hpp>
#include <windows.h>

using namespace std;

RESULT_CODE getMarkupData(SMarkup* cMarkupData, FILE* markupFile)
{
    int i;
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


void writeProjectionFeatureEntry(FILE* featureFile, SMarkup* cMarkup, int* pProjR, int* pProjL, int maxRadius, SPupilInfo psPI, SIrisInfo psII)
{
    // writing file name, CenX, CenY and projection length
    fprintf(featureFile, "%s,%d,%d,%d,", cMarkup->filename, cMarkup->pupil.xc, cMarkup->pupil.yc, maxRadius);
    //writing features
    for (int i = 0; i < maxRadius; ++i)
        fprintf(featureFile, "%d,", pProjR[i]);
    for (int i = 0; i < maxRadius; ++i)
        fprintf(featureFile, "%d,", pProjL[i]);
    // writing targets
    fprintf(featureFile, "%d,%d,%d,%d\n", cMarkup->BazRad.pu_r, cMarkup->BazRad.ir_r, psPI.r, psII.r);
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
    while ((res = getMarkupData(&cMarkupData, markupFile)) == ERROR_OK)
    {
        memset(imagePath, 0, MAX_FILENAME * sizeof(char));
        strcpy(imagePath, srcFolder);

        strcat(imagePath, cMarkupData.filename);
        printf("image path %s\n", imagePath);
        //Read the image data
        cv::Mat srcImage = cv::imread(imagePath, CV_8U);
        int H = srcImage.rows, W = srcImage.cols;
        printf("Markup entry %d: %s\n", i, cMarkupData.filename);
        printf("CenX: %d CenY: %d\n", cMarkupData.pupil.xc, cMarkupData.pupil.yc);

        int *pProjR = (int*)malloc(W * sizeof(int));
        int *pProjL = (int*)malloc(W * sizeof(int));

        memset(pProjR, 0, W * sizeof(int));
        memset(pProjL, 0, W * sizeof(int));
        SPupilInfo psPI;
        psPI.r = -1;
        psPI.xc = -1;
        psPI.yc = -1;
        SIrisInfo psII;
        psII.r = -1;
        psII.xc = -1;
        psII.yc = -1;
        SCenterInfo psCI;
        psCI.xc = cMarkupData.pupil.xc;
        psCI.yc = cMarkupData.pupil.yc;

        int kernel3x3[] = { -3, 0, 3, -10, 0, 10, -3, 0, 3 };
        int maxR = MIA_min(MIA_max(MIA_min(psCI.xc, W - psCI.xc), MIA_min(psCI.yc, H - psCI.yc)), H / 2);
        std::cout << "Max R: " << maxR << std::endl;

        if ((res = myDetectBaseRadii(&psPI, &psII, &psCI, srcImage.data, W, H, kernel3x3, -30, 60,
            MINIMAL_PUPIL_RADIUS, maxR, 3, 0.95, pProjR, pProjL)))
        {
            fprintf(stderr, "[ERROR]: Projection radii detection failed.\n");
            return res;
        }

        /*if ((res = FindHoughProjection(pProjR, pProjL, -30, 60, kernel3x3, 4, 0.95, MINIMAL_PUPIL_RADIUS, maxR, srcImage.data, W, H, psCI.xc, psCI.yc)))
        {
            fprintf(stderr, "[ERROR]: Projection feature generation failed.\n");
            return res;
        }
        memset(pProjL, 0, MINIMAL_PUPIL_RADIUS * sizeof(int));
        memset(pProjR, 0, MINIMAL_PUPIL_RADIUS * sizeof(int));

        int mode = IVIR_ARRPOXIR_PRPRPR | IVIR_ARRPOXIR_USEPROJ; //  IVIR_ARRPOXIR_LEFEYE | IVIR_ARRPOXIR_PRPRPR;
        if ((res = IVIR_PrPu(
            &psPI, &psII, &psCI, srcImage.data, W, H,
            mode, -1, -1, pProjR, pProjL, NULL)))
        {
            fprintf(stderr, "[ERROR]: base radii estimations failed.\n");
            return res;
        }*/
        printf("Found radii: %d %d, markup: %d %d\n", psPI.r, psII.r, cMarkupData.BazRad.pu_r, cMarkupData.BazRad.ir_r);

        writeProjectionFeatureEntry(featureFile, &cMarkupData, pProjR, pProjL, W / 2, psPI, psII);

        //         if ((res = DetectBaseRadii(&psPI, &psII, (const SCenterInfo*)&psCI, srcImage.data, W, H, pProjR, pProjL)))
        //         {
        //             free(pProjL);
        //             free(pProjR);
        //             fprintf(stderr, "ERROR: GenerateProjectionFeatures() returned %d.\n", res);
        //             return res;
        //         }

        // writeProjectionFeatureEntry(featureFile, &cMarkupData, pProjR, pProjL, W/2, psPI, psII);
        // 
        //         cv::namedWindow("MarkupImage", CV_WINDOW_AUTOSIZE); //create a window with the name "MyWindow"
        //         cv::imshow("MyWindow", srcImage); //display the image which is stored in the 'img' in the "MyWindow" window
        // 
        //         cv::waitKey(0); //wait infinite time for a keypress
        //         cv::destroyWindow("MarkupImage"); //destroy the window with the name, "MyWindow"
        printf("\n");
        i++;
        free(pProjL);
        free(pProjR);
        srcImage.release();
        // Sleep(1000);
    }
    free(imagePath);
    fclose(featureFile);
    fclose(markupFile);
    printf("processing result: %d\n", res);

    return ERROR_OK;
}



int main()
{
    RESULT_CODE res = GenerateProjectionFeatures(
        "D:\\data\\mixed\\",
        "D:\\data\\mixed\\params_res.txt",
        "D:\\data\\BaseRadii\\ProjectionFeatures7WithRad.csv");

    printf("GenerateProjectionFeatures() returns %d.\n", res);

//     char str;
//     cin >> str;
	return 0;
}