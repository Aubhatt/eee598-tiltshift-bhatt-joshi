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

//#define WEIGHT_MATRIX
#define WEIGHT_VECTOR

void gKernel_create(float sigma, int twiceCeil_r,float *ptr_array)
{
    float twoSigmaSquare = 2*sigma*sigma;
    float sum=0;
    int r = twiceCeil_r/2;
    int x,y;
    float num;

    //LOGI("gKernel_create called \n");

    for(int i=0; i<=twiceCeil_r; i++)
    {
        for(int j=0; j<=twiceCeil_r; j++)
        {
            x = j-r;
            y = i-r;
            num = sqrt(x*x + y*y);
            ptr_array[i*twiceCeil_r+j] = exp(-(num*num)/twoSigmaSquare)/(M_PI*twoSigmaSquare);
            sum = sum + ptr_array[i*twiceCeil_r+j];
        }
    }

    for(int i=0; i<=twiceCeil_r; i++)
    {
        for(int j=0; j<=twiceCeil_r; j++)
        {
            ptr_array[i*twiceCeil_r+j] = ptr_array[i*twiceCeil_r+j]/sum;
            // LOGI("ptr_array[%d*%d+%d] = %f",i,twiceCeil_r,j,ptr_array[i*twiceCeil_r+j]);
        }
    }

}






void gKernel_create_vector(float sigma, int twiceCeil_r,float *ptr_array)
{
    float twoSigmaSquare = 2*sigma*sigma;
    float sum=0;
    int r = twiceCeil_r/2;
    int k;
    float num,den;
    LOGI("gKernel_create called \n");

    for(int j=0; j<=twiceCeil_r; j++)
    {
        k = j-r;
        num = sqrt(k*k);
        den = sqrt(twoSigmaSquare*M_PI);
        ptr_array[j] = exp(-(num*num)/twoSigmaSquare)/(den);
        sum = sum + ptr_array[j];
    }

    for(int j=0; j<=twiceCeil_r; j++) {
        ptr_array[j] = ptr_array[j] / sum;
        //LOGI("ptr_array[%d] = %f", j, ptr_array[j]);
    }

}




extern "C"
JNIEXPORT jint JNICALL



Java_edu_asu_ame_meteor_speedytiltshift2018_SpeedyTiltShift_tiltshiftcppnative(JNIEnv *env, jobject instance, jintArray inputPixels_, jintArray outputPixels_, jint width, jint height, jfloat sigma_far, jfloat sigma_near, jint a0, jint a1, jint a2, jint a3) {
    jint *pixels = env->GetIntArrayElements(inputPixels_, NULL);
    jint *outputPixels = env->GetIntArrayElements(outputPixels_, NULL);
    // LOGI("SpeedyTiltShift C++ Called \n");
0

#ifdef WEIGHT_VECTOR
    int ceil_rFar = ceil(2 * sigma_far);
    int ceil_rNear = ceil(2 * sigma_near);
    ceil_rNear = ceil_rFar = fmax(ceil_rFar,ceil_rNear);

    float gKernel_far[(2 * ceil_rFar + 1) * (2 * ceil_rFar + 1)];
    float gKernel_near[(2 * ceil_rNear + 1) * (2 * ceil_rNear + 1)];


    if (sigma_far > 0.6) {
        gKernel_create_vector(sigma_far, 2 * ceil_rFar, gKernel_far);
    } else {
        LOGI("SIGMA_FAR < 0.6");
    }

    if (sigma_near > 0.6) {
        gKernel_create_vector(sigma_near, 2 * ceil_rNear, gKernel_near);
    } else {
        LOGI("SIGMA_NEAR < 0.6");
    }

    int A = 0xFF;
    uint32_t R,G,B;
    uint32_t color;
    int r = ceil_rFar;
    float sum_r,sum_g,sum_b=0;
    int q[width*(height)];
    for(int j = ceil_rFar; j < a0; j++)
    {
        for (int i = ceil_rFar; i < width-ceil_rFar; i++)
        {
            for(int y = 0; y <= 2*ceil_rFar; y++)
            {
                R = (pixels[(j-r+y)*width+i] >> 16) & 0xFF;
                G = (pixels[(j-r+y)*width+i] >> 8) & 0xFF;
                B = (pixels[(j-r+y)*width+i]) & 0xFF;

                sum_b = sum_b + gKernel_far[y] * B;
                sum_g = sum_g + gKernel_far[y] * G;
                sum_r = sum_r + gKernel_far[y] * R;
            }
            B = (uint32_t) sum_b;
            G = (uint32_t) sum_g;
            R = (uint32_t) sum_r;
            color = (((A & 0xFF) << 24) | ((R & 0xFF) << 16) | ((G & 0xFF) << 8) | (B & 0xFF));
            q[j*width+i] = color;
            sum_b = sum_g = sum_r = 0;
        }
    }

    for(int j = ceil_rFar; j < a0; j++)
    {
        for (int i = ceil_rFar; i < width-ceil_rFar; i++)
        {
            for(int y = 0; y <= 2*ceil_rFar; y++)
            {
                R = (q[j*width+i-r+y] >> 16) & 0xFF;
                G = (q[j*width+i-r+y] >> 8) & 0xFF;
                B = q[j*width+i-r+y] & 0xFF;

                sum_b = sum_b + gKernel_far[y] * B;
                sum_g = sum_g + gKernel_far[y] * G;
                sum_r = sum_r + gKernel_far[y] * R;
            }
            B = (uint32_t) sum_b;
            G = (uint32_t) sum_g;
            R = (uint32_t) sum_r;
            color = (((A & 0xFF) << 24) | ((R & 0xFF) << 16) | ((G & 0xFF) << 8) | (B & 0xFF));
            outputPixels[j*width+i] = color;
            sum_b = sum_g = sum_r = 0;
        }
    }


    float sigma_dynamic = 1;
    float expres = 1;
    int numerator, denominator;
    float gKernel_far_two[(2 * ceil_rFar + 1) * (2 * ceil_rFar + 1)];
    for(int j = a0; j < a1; j++)
    {
        numerator = a1 - j;
        denominator = a1 - a0;
        expres = (float)numerator/denominator;
        sigma_dynamic = sigma_far * expres;
        gKernel_create_vector(sigma_dynamic,2*ceil_rFar,gKernel_far_two);
        for (int i = ceil_rFar; i < width - ceil_rFar; i++)
        {
            for(int y = 0; y <= 2*ceil_rFar; y++)
            {
                R = (pixels[(j-r+y)*width+i] >> 16) & 0xFF;
                G = (pixels[(j-r+y)*width+i] >> 8) & 0xFF;
                B = (pixels[(j-r+y)*width+i]) & 0xFF;


                sum_b = sum_b + gKernel_far_two[y] * B;
                sum_g = sum_g + gKernel_far_two[y] * G;
                sum_r = sum_r + gKernel_far_two[y] * R;
            }
            B = (uint32_t) sum_b;
            G = (uint32_t) sum_g;
            R = (uint32_t) sum_r;
            color = (((A & 0xFF) << 24) | ((R & 0xFF) << 16) | ((G & 0xFF) << 8) | (B & 0xFF));
            q[j*width+i] = color;
            sum_b = sum_g = sum_r = 0;
        }
    }

    for(int j = a0; j < a1; j++)
    {
        for (int i = ceil_rFar; i < width-ceil_rFar; i++)
        {
            for(int y = 0; y <= 2*ceil_rFar; y++)
            {
                R = (q[j*width+i-r+y] >> 16) & 0xFF;
                G = (q[j*width+i-r+y] >> 8) & 0xFF;
                B = q[j*width+i-r+y] & 0xFF;

                sum_b = sum_b + gKernel_far_two[y] * B;
                sum_g = sum_g + gKernel_far_two[y] * G;
                sum_r = sum_r + gKernel_far_two[y] * R;
            }
            B = (uint32_t) sum_b;
            G = (uint32_t) sum_g;
            R = (uint32_t) sum_r;
            color = (((A & 0xFF) << 24) | ((R & 0xFF) << 16) | ((G & 0xFF) << 8) | (B & 0xFF));
            outputPixels[j*width+i] = color;
            sum_b = sum_g = sum_r = 0;
        }
    }

    for(int j = a1; j < a2; j++)
    {
        for (int i = ceil_rFar; i < width - ceil_rFar ;i++)
        {
            outputPixels[j*width + i] = pixels[j*width + i];
        }
    }

    float gKernel_near_two[(2 * ceil_rNear + 1) * (2 * ceil_rNear + 1)];
    for(int j=a2; j < a3; j++)
    {
        numerator = j - a2;
        denominator = a3 - a2;
        expres = (float)numerator/denominator;
        sigma_dynamic = sigma_near * expres;
        gKernel_create_vector(sigma_dynamic,2*ceil_rNear,gKernel_near_two);
        for (int i = ceil_rNear; i < width - ceil_rNear; i++)
        {
            for(int y = 0; y <= 2*ceil_rFar; y++)
            {
                R = (pixels[(j-r+y)*width+i] >> 16) & 0xFF;
                G = (pixels[(j-r+y)*width+i] >> 8) & 0xFF;
                B = (pixels[(j-r+y)*width+i]) & 0xFF;


                sum_b = sum_b + gKernel_near_two[y] * B;
                sum_g = sum_g + gKernel_near_two[y] * G;
                sum_r = sum_r + gKernel_near_two[y] * R;
            }
            B = (uint32_t) sum_b;
            G = (uint32_t) sum_g;
            R = (uint32_t) sum_r;
            color = (((A & 0xFF) << 24) | ((R & 0xFF) << 16) | ((G & 0xFF) << 8) | (B & 0xFF));
            q[j*width+i] = color;
            sum_b = sum_g = sum_r = 0;
        }
    }

    for(int j = a2; j < a3; j++)
    {
        for (int i = ceil_rFar; i < width - ceil_rNear; i++)
        {
            for (int y = 0; y <= 2 * ceil_rNear; y++)
            {
                R = (q[j*width+i-r+y] >> 16) & 0xFF;
                G = (q[j*width+i-r+y] >> 8) & 0xFF;
                B = q[j*width+i-r+y] & 0xFF;

                sum_b = sum_b + gKernel_near_two[y] * B;
                sum_g = sum_g + gKernel_near_two[y] * G;
                sum_r = sum_r + gKernel_near_two[y] * R;
            }
            B = (uint32_t) sum_b;
            G = (uint32_t) sum_g;
            R = (uint32_t) sum_r;
            color = (((A & 0xFF) << 24) | ((R & 0xFF) << 16) | ((G & 0xFF) << 8) | (B & 0xFF));
            outputPixels[j*width+i] = color;
            sum_b = sum_g = sum_r = 0;
        }
    }


    for(int j = a3; j < height-ceil_rNear; j++)
    {
        for (int i = ceil_rNear; i < width - ceil_rNear; i++)
        {
            for(int y = 0; y <= 2*ceil_rNear; y++)
            {
                R = (pixels[(j-r+y)*width+i] >> 16) & 0xFF;
                G = (pixels[(j-r+y)*width+i] >> 8) & 0xFF;
                B = (pixels[(j-r+y)*width+i]) & 0xFF;

                sum_b = sum_b + gKernel_near[y] * B;
                sum_g = sum_g + gKernel_near[y] * G;
                sum_r = sum_r + gKernel_near[y] * R;
            }
            B = (uint32_t) sum_b;
            G = (uint32_t) sum_g;
            R = (uint32_t) sum_r;
            color = (((A & 0xFF) << 24) | ((R & 0xFF) << 16) | ((G & 0xFF) << 8) | (B & 0xFF));
            q[j*width+i] = color;
            sum_b = sum_g = sum_r = 0;
        }
    }

    for(int j = a3; j < height-ceil_rNear; j++)
    {
        for (int i = ceil_rNear; i < width - ceil_rNear; i++)
        {
            for (int y = 0; y <= 2 * ceil_rFar; y++) {
                R = (q[j*width+i-r+y] >> 16) & 0xFF;
                G = (q[j*width+i-r+y] >> 8) & 0xFF;
                B = q[j*width+i-r+y] & 0xFF;

                sum_b = sum_b + gKernel_near[y] * B;
                sum_g = sum_g + gKernel_near[y] * G;
                sum_r = sum_r + gKernel_near[y] * R;
            }
            B = (uint32_t) sum_b;
            G = (uint32_t) sum_g;
            R = (uint32_t) sum_r;
            color = (((A & 0xFF) << 24) | ((R & 0xFF) << 16) | ((G & 0xFF) << 8) | (B & 0xFF));
            outputPixels[j * width + i] = color;
            sum_b = sum_g = sum_r = 0;
        }
    }
#endif



#ifdef WEIGHT_MATRIX
    int ceil_rFar = ceil(2 * sigma_far);
    int ceil_rNear = ceil(2 * sigma_near);
    ceil_rNear = ceil_rFar = fmax(ceil_rFar,ceil_rNear);


    float gKernel_far[(2 * ceil_rFar + 1) * (2 * ceil_rFar + 1)];
    float gKernel_near[(2 * ceil_rNear + 1) * (2 * ceil_rNear + 1)];



    if (sigma_far > 0.6) {
        gKernel_create(sigma_far, 2*ceil_rFar, gKernel_far);
    }
    if (sigma_near > 0.6) {
        gKernel_create(sigma_near, 2*ceil_rNear, gKernel_near);
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

    float sigma_dynamic = 1;
    float expres = 1;
    int numerator, denominator;
    float gKernel_far_two[(2 * ceil_rFar + 1) * (2 * ceil_rFar + 1)];
    for (int j = a0; j < a1; j++)    // Row of original image
    {
        numerator = a1 - j;
        denominator = a1 - a0;
        expres = (float)numerator/denominator;
        sigma_dynamic = sigma_far * expres;
        gKernel_create(sigma_dynamic,2*ceil_rFar,gKernel_far_two);
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
        for (int i = ceil_rFar; i < width - ceil_rFar ;i++)
        {
            outputPixels[j*width + i] = pixels[j*width + i];
        }
    }


    float gKernel_near_two[(2 * ceil_rNear + 1) * (2 * ceil_rNear + 1)];
    for (int j=a2; j < a3; j++)    // Row of original image
    {
        numerator = j - a2;
        denominator = a3 - a2;
        expres = (float)numerator/denominator;
        sigma_dynamic = sigma_far * expres;
        gKernel_create(sigma_dynamic,2*ceil_rFar,gKernel_near_two);
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
            for (int yG = 0; yG <= 2 * ceil_rNear; yG++)  // Row of Kernel
            {
                for (int xG = 0; xG <= 2 * ceil_rNear; xG++)  // Column of Kernel
                {
                    x = i - ceil_rNear + xG;
                    y = j - ceil_rNear + yG;
                    B = (pixels[y * width + x]) & 0xFF;
                    G = (pixels[y * width + x] >> 8) & 0xFF;
                    R = (pixels[y * width + x] >> 16) & 0xFF;
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
            outputPixels[j * width + i] = color;
            sum_b = sum_g = sum_r = 0;

        }
    }

#endif

    //LOGI("DONE \n");
    env->ReleaseIntArrayElements(inputPixels_, pixels, 0);
    env->ReleaseIntArrayElements(outputPixels_, outputPixels, 0);
    return 0;
}
