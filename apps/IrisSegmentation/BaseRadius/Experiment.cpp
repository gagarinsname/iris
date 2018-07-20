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
            MINIMAL_PUPIL_RADIUS, maxR, 3, 0.95, pProjR, pProjL, 4)))
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
    /*RESULT_CODE res = GenerateProjectionFeatures(
        "D:\\data\\mixed\\",
        "D:\\data\\mixed\\params_res.txt",
        "D:\\data\\BaseRadii\\ProjectionFeatures7WithRad.csv");*/
	char imagePath[] = "D:\\data\\IrisDB\\casia3-lamp\\2244(R-1-1)060805-082000_e_01_0002.bmp";
	cv::Mat srcImage = cv::imread(imagePath, CV_8U);
	int H = srcImage.rows, W = srcImage.cols;
	// printf("Markup entry %d: %s\n", i, cMarkupData.filename);
	// printf("CenX: %d CenY: %d\n", cMarkupData.pupil.xc, cMarkupData.pupil.yc);

	int *pProjR = (int*)malloc(W * sizeof(int));
	int *pProjL = (int*)malloc(W * sizeof(int));

	memset(pProjR, 0, W * sizeof(int));
	memset(pProjL, 0, W * sizeof(int));

	// int kernel3x3[] = { -1, 0, 1, -2, 0, 2, -1, 0, 1 };
	// int res = FindHoughProjection(pProjR, pProjL, -30, 50, kernel3x3, 4, 0.95, 10, 320, srcImage.data, W, H, 356, 222);

    // printf("GenerateProjectionFeatures() returns %d.\n", res);
//     char str;
//     cin >> str;
	int a[] = { 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		0,    0,    0,    0,    0, 1064, 1429, 1813, 2222, 2273, 2419,
		2393, 2078, 1737, 1789, 1688, 1552, 1370, 1601, 1676, 1747, 1886,
		2122, 2147, 2073, 2158, 2357, 2481, 2576, 2615, 2811, 2953, 3008,
		3159, 3154, 2962, 2981, 3227, 3470, 3684, 3748, 3763, 3646, 3464,
		3442, 3153, 2685, 2279, 1856, 1505, 1225, 1102, 1102, 1043, 1079,
		1126, 1059,  985,  963,  987, 1002,  921,  886,  830,  795,  831,
		869,  893,  926,  899,  891,  802,  700,  704,  693,  721,  664,
		571,  542,  580,  644,  704,  673,  733,  801,  905,  919,  988,
		1014,  982, 1004, 1004,  940,  809,  829,  905,  900,  854,  885,
		976, 1090, 1154, 1238, 1198, 1242, 1293, 1388, 1473, 1368, 1270,
		1344, 1465, 1689, 1975, 2330, 2725, 3096, 3569, 4062, 4418, 4729,
		4904, 4974, 4944, 4804, 4584, 4349, 4108, 3788, 3444, 3070, 2650,
		2279, 1959, 1717, 1535, 1351, 1299, 1289, 1291, 1349, 1384, 1386,
		1388, 1400, 1398, 1324, 1191, 1187, 1196, 1238, 1340, 1467, 1496,
		1516, 1536, 1627, 1681, 1581, 1534, 1444, 1337, 1280, 1209, 1205,
		1176, 1170, 1272, 1214, 1186, 1185, 1201, 1228, 1142, 1165, 1127,
		1038, 1007,  991,  970,  920,  909,  922,  858,  740,  808,  891,
		969, 1036 };

	int b[] = { 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
		0,    0,    0,    0,    0,  449,  385,  280,  178,  177,  253,
		312,  382,  417,  526,  607,  725,  869,  989, 1012, 1035, 1084,
		1141, 1146, 1169, 1189, 1198, 1130, 1192, 1120, 1029, 1029, 1097,
		1120, 1112, 1077, 1065, 1091, 1396, 1798, 2132, 2350, 2450, 2478,
		2431, 2445, 2342, 2056, 1643, 1287, 1017,  920,  918,  993, 1126,
		1202, 1235, 1269, 1296, 1341, 1349, 1305, 1276, 1106, 1012,  942,
		884,  798,  694,  614,  597,  561,  579,  633,  714,  830,  886,
		912,  931,  889,  911,  989, 1012, 1032,  970,  981,  997,  976,
		987,  967,  925,  894,  903,  935,  943,  940, 1006, 1064, 1106,
		1175, 1267, 1380, 1626, 1949, 2374, 2790, 3243, 3685, 4100, 4489,
		4789, 4917, 4983, 4943, 4894, 4779, 4669, 4510, 4338, 4186, 4056,
		3901, 3750, 3568, 3416, 3220, 3069, 2925, 2733, 2627, 2540, 2426,
		2287, 2144, 2056, 1969, 1898, 1837, 1706, 1570, 1454, 1431, 1416,
		1355, 1299, 1210, 1127, 1060, 1029, 1031, 1003,  936,  923,  858,
		834,  868,  915,  886,  849,  798,  801,  774,  869,  917,  886,
		887,  961,  977,  992, 1034, 1072, 1041, 1044, 1123, 1135, 1132,
		1179, 1215, 1216, 1258, 1263, 1241, 1206, 1190, 1224, 1254, 1309,
		1352, 1372 };
	
	int sizeR = 200;
	int sizeP = 200 * 2 + 2;
	sPnt* path = (sPnt*)malloc(sizeP * sizeof(sPnt));
	float *matDTW = (float*)malloc(sizeR * sizeR * sizeof(float));
	int nScales = 2;
	// int res = myDynamicTimeWarpingPyramid(path, &sizeP, matDTW, a, 38, b, 38, 3, 1, &nScales);
	int res = calcFullDynamicTimeWarpingWrapper(path, &sizeP, matDTW, b, 200, a, 200, -1, 2.0f);
	printf("calcFullDynamicTimeWarpingWrapper() returns %d.\n", res);

	free(matDTW);
	free(path);
	free(pProjR);
	free(pProjL);
	return 0;
}