#include <jni.h>
#include <string>
#include <cpu-features.h>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <android/log.h>

#define TAG "SPEEDY_TS"
#define LOGE(...) __android_log_print(ANDROID_LOG_ERROR, TAG, __VA_ARGS__)
#define LOGW(...) __android_log_print(ANDROID_LOG_WARN, TAG, __VA_ARGS__)
#define LOGI(...) __android_log_print(ANDROID_LOG_INFO, TAG, __VA_ARGS__)
#define LOGD(...) __android_log_print(ANDROID_LOG_DEBUG, TAG, __VA_ARGS__)


void gKernel_create(float sigma, int twiceCeil_r,float *ptr_array)
{
    float twoSigmaSquare = 2*sigma*sigma;
    float sum=0;
    int r = twiceCeil_r/2;
    int x,y;
    float num;
    LOGI("gKernel_create called \n");

    for(int i=0; i<=twiceCeil_r; i++)
    {
        for(int j=0; j<=twiceCeil_r; j++)
        {
            x = j-r;
            y = i-r;
            num = sqrt(x*x + y*y);
            std::cout << exp(-(num*num)/twoSigmaSquare)/(M_PI*twoSigmaSquare) << "\n";
            ptr_array[i*twiceCeil_r+j] = exp(-(num*num)/twoSigmaSquare)/(M_PI*twoSigmaSquare);
            sum = sum + ptr_array[i*twiceCeil_r+j];
        }
    }

    for(int i=0; i<=twiceCeil_r; i++)
    {
        for(int j=0; j<=twiceCeil_r; j++)
        {
            ptr_array[i*twiceCeil_r+j] = ptr_array[i*twiceCeil_r+j]/sum;

        }
    }

}

void gKernel_create_two_far(float sigma, int twiceCeil_r,float *ptr_array,int a0,int a1)
{
    float twoSigmaSquare;
    float func_sigma;
    float sum=0;
    int r = twiceCeil_r/2;
    int x,y;
    int y_sigma;
    float var,num;
    LOGI("gKernel_create called \n");

    for(int i=0; i<=twiceCeil_r; i++)
    {
        for(int j=0; j<=twiceCeil_r; j++)
        {
            x = j-r;
            y = i-r;
            num = sqrt(x*x + y*y);
            y_sigma = i;
            var = (a1-y_sigma)/(a1-a0);
            func_sigma = sigma * var;
            twoSigmaSquare = 2 * func_sigma * func_sigma;
            ptr_array[i*twiceCeil_r+j] = exp(-(num*num)/twoSigmaSquare)/(M_PI*twoSigmaSquare);
            sum = sum + ptr_array[i*twiceCeil_r+j];
        }

    }

    for(int i=0; i<=twiceCeil_r; i++)
    {
        for(int j=0; j<=twiceCeil_r; j++)
        {
            ptr_array[i*twiceCeil_r+j] = ptr_array[i*twiceCeil_r+j]/sum;

        }
    }

}


void gKernel_create_two_near(float sigma, int twiceCeil_r,float *ptr_array,int a2,int a3)
{
    float twoSigmaSquare;
    float func_sigma;
    float sum=0;
    int r = twiceCeil_r/2;
    int x,y;
    int y_sigma;
    float var,num;
    LOGI("gKernel_create called \n");

    for(int i=0; i<=twiceCeil_r; i++)
    {
        for(int j=0; j<=twiceCeil_r; j++)
        {
            x = j-r;
            y = i-r;
            num = sqrt(x*x + y*y);
            y_sigma = i;
            var = (y_sigma-a2)/(a3-a2);
            func_sigma = sigma * var;
            twoSigmaSquare = 2 * func_sigma * func_sigma;
            ptr_array[i*twiceCeil_r+j] = exp(-(num*num)/twoSigmaSquare)/(M_PI*twoSigmaSquare);
            sum = sum + ptr_array[i*twiceCeil_r+j];
        }

    }

    for(int i=0; i<=twiceCeil_r; i++)
    {
        for(int j=0; j<=twiceCeil_r; j++)
        {
            ptr_array[i*twiceCeil_r+j] = ptr_array[i*twiceCeil_r+j]/sum;
            LOGI("ptr_array[%d*%d+%d] = %f",i,twiceCeil_r,j,ptr_array[i*twiceCeil_r+j]);
        }
    }

}



extern "C"
JNIEXPORT jint JNICALL



Java_edu_asu_ame_meteor_speedytiltshift2018_SpeedyTiltShift_tiltshiftcppnative(JNIEnv *env, jobject instance, jintArray inputPixels_, jintArray outputPixels_, jint width, jint height, jfloat sigma_far, jfloat sigma_near, jint a0, jint a1, jint a2, jint a3) {
    jint *pixels = env->GetIntArrayElements(inputPixels_, NULL);
    jint *outputPixels = env->GetIntArrayElements(outputPixels_, NULL);
    LOGI("SpeedyTiltShift C++ Called \n");



    int ceil_rFar = ceil(2 * sigma_far);
    int ceil_rNear = ceil(2 * sigma_near);


    float gKernel_far[(2 * ceil_rFar + 1) * (2 * ceil_rFar + 1)];
    float gKernel_near[(2 * ceil_rNear + 1) * (2 * ceil_rNear + 1)];

    float gKernel_far_two[(2 * ceil_rFar + 1) * (2 * ceil_rFar + 1)];
    float gKernel_near_two[(2 * ceil_rNear + 1) * (2 * ceil_rNear + 1)];

    if (sigma_far > 0.6) {
        gKernel_create(sigma_far, 2*ceil_rFar, gKernel_far);
        gKernel_create_two_far(sigma_far,2*ceil_rFar,gKernel_far_two,a0,a1);
    }
    if (sigma_near > 0.6) {
        gKernel_create(sigma_near, 2*ceil_rNear, gKernel_near);
        gKernel_create_two_near(sigma_near,2*ceil_rNear,gKernel_near_two,a2,a3);
    }

    float sum_r, sum_g, sum_b = 0;
    float value_r, value_g, value_b = 0;
    uint32_t color;
    uint32_t A = 0xFF;
    uint32_t  R, G, B;
    int x,y;


    for (int j = ceil_rFar; j < a0; j++)    // Row of original image
    {
        for (int i = ceil_rFar; i < width - ceil_rFar; i++)  // Column of original image
        {
            for (int yG = 0; yG <= 2*ceil_rFar; yG++)  // Row of Kernel
            {
                for (int xG = 0; xG <= 2*ceil_rFar; xG++)  // Column of Kernel
                {
                    x = i - ceil_rFar + xG;
                    y = j - ceil_rFar + yG;
                    B = (pixels[y*width+x]) & 0xFF;
                    G = (pixels[y*width+x]>>8) & 0xFF;
                    R = (pixels[y*width+x]>>16) & 0xFF;
                    value_b = gKernel_far[yG * 2 * ceil_rFar + xG] * B;
                    sum_b = sum_b + value_b;

                    value_g = gKernel_far[yG * 2 * ceil_rFar + xG] * G;
                    sum_g = sum_g + value_g;

                    value_r = gKernel_far[yG * 2 * ceil_rFar + xG] * R;
                    sum_r = sum_r + value_r;
                }

            }
            B = (uint32_t) sum_b;
            G = (uint32_t) sum_g;
            R = (uint32_t) sum_r;
            color = (((A & 0xFF) << 24) | ((R & 0xFF) << 16) | ((G & 0xFF) << 8) | (B & 0xFF));
            outputPixels[j*width + i] = color;
            sum_b = sum_g = sum_r = 0;

        }
    }


    for (int j = a0; j < a1; j++)    // Row of original image
    {
        for (int i = ceil_rFar; i < width - ceil_rFar; i++)  // Column of original image
        {
            for (int yG = 0; yG <= 2*ceil_rFar; yG++)  // Row of Kernel
            {
                for (int xG = 0; xG <= 2*ceil_rFar; xG++)  // Column of Kernel
                {
                    x = i - ceil_rFar + xG;
                    y = j - ceil_rFar + yG;
                    B = (pixels[y*width+x]) & 0xFF;
                    G = (pixels[y*width+x]>>8) & 0xFF;
                    R = (pixels[y*width+x]>>16) & 0xFF;
                    value_b = gKernel_far_two[yG * 2 * ceil_rFar + xG] * B;
                    sum_b = sum_b + value_b;

                    value_g = gKernel_far_two[yG * 2 * ceil_rFar + xG] * G;
                    sum_g = sum_g + value_g;

                    value_r = gKernel_far_two[yG * 2 * ceil_rFar + xG] * R;
                    sum_r = sum_r + value_r;
                }
            }
            B = (uint32_t) sum_b;
            G = (uint32_t) sum_g;
            R = (uint32_t) sum_r;
            color = (((A & 0xFF) << 24) | ((R & 0xFF) << 16) | ((G & 0xFF) << 8) | (B & 0xFF));
            outputPixels[j*width + i] = color;
            sum_b = sum_g = sum_r = 0;

        }
    }

    for(int j = a1; j < a2; j++)
    {
        for (int i = 0; i < width;i++)
        {
            outputPixels[j*width + i] = pixels[j*width + i];
        }
    }



    for (int j=a2; j < a3; j++)    // Row of original image
    {
        for (int i = ceil_rNear; i < width - ceil_rNear; i++)  // Column of original image
        {
            for (int yG = 0; yG <= 2*ceil_rNear; yG++)  // Row of Kernel
            {
                for (int xG = 0; xG <= 2*ceil_rNear; xG++)  // Column of Kernel
                {

                    x = i - ceil_rNear + xG;
                    y = j - ceil_rNear + yG;
                    B = (pixels[y*width+x]) & 0xFF;
                    G = (pixels[y*width+x]>>8) & 0xFF;
                    R = (pixels[y*width+x]>>16) & 0xFF;
                    value_b = gKernel_near_two[yG * 2 * ceil_rNear + xG] * B;
                    sum_b = sum_b + value_b;

                    value_g = gKernel_near_two[yG * 2 * ceil_rNear + xG] * G;
                    sum_g = sum_g + value_g;

                    value_r = gKernel_near_two[yG * 2 * ceil_rNear + xG] * R;
                    sum_r = sum_r + value_r;
                }

            }
            B = (uint32_t) sum_b;
            G = (uint32_t) sum_g;
            R = (uint32_t) sum_r;
            color = (((A & 0xFF) << 24) | ((R & 0xFF) << 16) | ((G & 0xFF) << 8) | (B & 0xFF));
            outputPixels[j*width + i] = color;
            sum_b = sum_g = sum_r = 0;

        }
    }



    for (int j=a3;j<height-ceil_rNear;j++)    // Row of original image
    {
        for (int i = ceil_rNear; i < width - ceil_rNear; i++)  // Column of original image
        {
            for (int yG = 0; yG <= 2*ceil_rNear; yG++)  // Row of Kernel
            {
                for (int xG = 0; xG <= 2*ceil_rNear; xG++)  // Column of Kernel
                {
                    x = i - ceil_rNear + xG;
                    y = j - ceil_rNear + yG;
                    B = (pixels[y*width+x]) & 0xFF;
                    G = (pixels[y*width+x]>>8) & 0xFF;
                    R = (pixels[y*width+x]>>16) & 0xFF;
                    value_b = gKernel_near[yG * 2 * ceil_rNear + xG] * B;
                    sum_b = sum_b + value_b;

                    value_g = gKernel_near[yG * 2 * ceil_rNear + xG] * G;
                    sum_g = sum_g + value_g;

                    value_r = gKernel_near[yG * 2 * ceil_rNear + xG] * R;
                    sum_r = sum_r + value_r;
                }

            }
            B = (uint32_t) sum_b;
            G = (uint32_t) sum_g;
            R = (uint32_t) sum_r;
            color = (((A & 0xFF) << 24) | ((R & 0xFF) << 16) | ((G & 0xFF) << 8) | (B & 0xFF));
            outputPixels[j*width + i] = color;
            sum_b = sum_g = sum_r = 0;

        }
    }

    LOGI("DONE \n");
    env->ReleaseIntArrayElements(inputPixels_, pixels, 0);
    env->ReleaseIntArrayElements(outputPixels_, outputPixels, 0);
    return 0;
}

