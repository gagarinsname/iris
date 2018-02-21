#include <stdio.h>
#include <string.h>
#include <opencv2/opencv.hpp>

#include "GlobalParams.h"
#include "BaseRadii.h"
using namespace std;

int main()
{
    RESULT_CODE res = GenerateProjectionFeatures(
        "D:\\data\\mixed\\",
        "D:\\data\\mixed\\params_res.txt",
        "D:\\data\\ProjectionFeatures6.csv");

    printf("GenerateProjectionFeatures() returns %d.\n", res);

    char str;
    cin << str;
	return 0;
}