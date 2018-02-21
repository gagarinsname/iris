#include <stdio.h>
#include <string.h>
#include <opencv2/opencv.hpp>

#include "GlobalParams.h"
#include "BaseRadii.h"

int main()
{
    RESULT_CODE res = GenerateProjectionFeatures(
        "D:\\data\\mixed\\",
        "D:\\data\\mixed\\params_res.txt",
        "D:\\data\\ProjectionFeatures6.csv");

    printf("GenerateProjectionFeatures() returned %d", res);

    std::string str;
    std::getline(std::cin, str);
	return 0;
}