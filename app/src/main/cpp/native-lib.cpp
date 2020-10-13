#include <jni.h>
#include <string>
#include <cpu-features.h>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <android/log.h>
#include <arm_neon.h>
//#include <armintr.h>
//#include <arm64intr.h>
#define TAG "SPEEDY_TS"
#define LOGE(...) __android_log_print(ANDROID_LOG_ERROR, TAG, __VA_ARGS__)
#define LOGW(...) __android_log_print(ANDROID_LOG_WARN, TAG, __VA_ARGS__)
#define LOGI(...) __android_log_print(ANDROID_LOG_INFO, TAG, __VA_ARGS__)
#define LOGD(...) __android_log_print(ANDROID_LOG_DEBUG, TAG, __VA_ARGS__)

//#define WEIGHT_MATRIX
//#define WEIGHT_VECTOR

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


extern "C"
JNIEXPORT jint JNICALL
Java_edu_asu_ame_meteor_speedytiltshift2018_SpeedyTiltShift_tiltshiftneonnative(JNIEnv *env,
                                                                                jclass instance,
                                                                                jintArray inputPixels_,
                                                                                jintArray outputPixels_,
                                                                                jint width,
                                                                                jint height,
                                                                                jfloat sigma_far,
                                                                                jfloat sigma_near,
                                                                                jint a0, jint a1,
                                                                                jint a2, jint a3) {
    jint *pixels = env->GetIntArrayElements(inputPixels_, NULL);
    jint *outputPixels = env->GetIntArrayElements(outputPixels_, NULL);


    uint8_t * arrayInPtr = (uint8_t *)pixels;
    uint8_t * arrayOutPtr = (uint8_t *)outputPixels;

    uint8_t q[width*height*4];
    uint8_t *arrayQptr = (uint8_t *)q;

    int ceil_rFar = ceil(2 * sigma_far);
    int ceil_rNear = ceil(2 * sigma_near);
    ceil_rNear = ceil_rFar = fmax(ceil_rFar,ceil_rNear);

    float gKernel_far[(2 * ceil_rFar + 1) ];
    float gKernel_near[(2 * ceil_rNear + 1)];
    float gKernel_far_two[(2 * ceil_rFar + 1)];
    float gKernel_near_two[(2 * ceil_rFar + 1) ];


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

    int var = height*width;
    for(int iter=0; iter<=2*ceil_rFar; iter++)
    {
        for (int i=0; i <width*height; i+= 16){
            uint8x16x4_t pixelChannels = vld4q_u8(arrayInPtr);
            int y = i/width;

            uint8x16_t Bvector = pixelChannels.val[0];
            uint8x16_t Gvector = pixelChannels.val[1];
            uint8x16_t Rvector = pixelChannels.val[2];
            uint8x16_t Avector = pixelChannels.val[3];

            uint8x8_t Blow = vget_low_u8(Bvector);
            uint8x8_t Glow = vget_low_u8(Gvector);
            uint8x8_t Rlow = vget_low_u8(Rvector);

            uint16x8_t Blow16 = vmovl_u8(Blow);
            uint16x8_t Glow16 = vmovl_u8(Glow);
            uint16x8_t Rlow16 = vmovl_u8(Rlow);

            if(y<a0)
            {
                if(sigma_far > 0.6)
                {
                    Blow16 = vmulq_n_u16(Blow16,(uint16_t)(gKernel_far[iter]*64));
                    Glow16 = vmulq_n_u16(Glow16,(uint16_t)(gKernel_far[iter]*64));
                    Rlow16 = vmulq_n_u16(Rlow16,(uint16_t)(gKernel_far[iter]*64));
                }
            }
            else if(y >= a0 && y < a1)
            {
                float sigmaFar_2 = (a1-y)/(a1-a0);
                sigmaFar_2 = sigma_far*sigmaFar_2;
                if(sigmaFar_2 > 0.6)
                {
                    gKernel_create_vector(sigmaFar_2,2*ceil_rFar,gKernel_far_two);
                }
                Blow16 = vmulq_n_u16(Blow16,(uint16_t)(gKernel_far[iter]*64));
                Glow16 = vmulq_n_u16(Glow16,(uint16_t)(gKernel_far[iter]*64));
                Rlow16 = vmulq_n_u16(Rlow16,(uint16_t)(gKernel_far[iter]*64));
            }
            else if(y >= a2 && y < a3)
            {
                float sigmaNear_2 = (y-a2)/(a2-a1);
                sigmaNear_2 = sigmaNear_2*sigma_near;
                if(sigmaNear_2 > 0.6)
                {
                    gKernel_create_vector(sigmaNear_2,2*ceil_rNear,gKernel_near_two);
                }
                Blow16 = vmulq_n_u16(Blow16,(uint16_t)(gKernel_near[iter]*64));
                Glow16 = vmulq_n_u16(Glow16,(uint16_t)(gKernel_near[iter]*64));
                Rlow16 = vmulq_n_u16(Rlow16,(uint16_t)(gKernel_near[iter]*64));
            }

            else
            {
                if(sigma_near > 0.6)
                {
                    Blow16 = vmulq_n_u16(Blow16,(uint16_t)(gKernel_near[iter]*64));
                    Glow16 = vmulq_n_u16(Glow16,(uint16_t)(gKernel_near[iter]*64));
                    Rlow16 = vmulq_n_u16(Rlow16,(uint16_t)(gKernel_near[iter]*64));
                }
            }




            Blow16 = vshrq_n_u16(Blow16,6);
            Glow16 = vshrq_n_u16(Rlow16,6);
            Rlow16 = vshrq_n_u16(Rlow16,6);

            Blow = vqmovn_u16(Blow16);
            Glow = vqmovn_u16(Glow16);
            Rlow = vqmovn_u16(Rlow16);

            uint8x8_t Bhigh = vget_high_u8(Bvector);
            uint8x8_t Ghigh = vget_high_u8(Gvector);
            uint8x8_t Rhigh = vget_high_u8(Rvector);

            uint16x8_t Bhigh16 = vmovl_u8(Bhigh);
            uint16x8_t Ghigh16 = vmovl_u8(Ghigh);
            uint16x8_t Rhigh16 = vmovl_u8(Rhigh);

            if(y<a0)
            {
                if(sigma_far > 0.6)
                {
                    Bhigh16 = vmulq_n_u16(Bhigh16,(uint16_t)(gKernel_far[iter]*64));
                    Ghigh16 = vmulq_n_u16(Ghigh16,(uint16_t)(gKernel_far[iter]*64));
                    Rhigh16 = vmulq_n_u16(Rhigh16,(uint16_t)(gKernel_far[iter]*64));
                }
            }
            else if(y >= a0 && y < a1)
            {
                float sigmaFar_2 = (a1-y)/(a1-a0);
                sigmaFar_2 = sigma_far*sigmaFar_2;
                if(sigmaFar_2 > 0.6)
                {
                    gKernel_create_vector(sigmaFar_2,2*ceil_rFar,gKernel_far_two);
                }
                Bhigh16 = vmulq_n_u16(Bhigh16,(uint16_t)(gKernel_far[iter]*64));
                Ghigh16 = vmulq_n_u16(Ghigh16,(uint16_t)(gKernel_far[iter]*64));
                Rhigh16 = vmulq_n_u16(Rhigh16,(uint16_t)(gKernel_far[iter]*64));
            }
            else if(y >= a2 && y < a3)
            {
                float sigmaNear_2 = (y-a2)/(a2-a1);
                sigmaNear_2 = sigmaNear_2*sigma_near;
                if(sigmaNear_2 > 0.6)
                {
                    gKernel_create_vector(sigmaNear_2,2*ceil_rNear,gKernel_near_two);
                }
                Bhigh16 = vmulq_n_u16(Bhigh16,(uint16_t)(gKernel_near[iter]*64));
                Ghigh16 = vmulq_n_u16(Ghigh16,(uint16_t)(gKernel_near[iter]*64));
                Rhigh16 = vmulq_n_u16(Rhigh16,(uint16_t)(gKernel_near[iter]*64));
            }
            else
            {
                if(sigma_near > 0.6)
                {
                    Bhigh16 = vmulq_n_u16(Bhigh16,(uint16_t)(gKernel_near[iter]*64));
                    Ghigh16 = vmulq_n_u16(Ghigh16,(uint16_t)(gKernel_near[iter]*64));
                    Rhigh16 = vmulq_n_u16(Rhigh16,(uint16_t)(gKernel_near[iter]*64));
                }

            }


            Bhigh16 = vmulq_n_u16(Bhigh16,(uint16_t)(gKernel_far[iter]*64));
            Ghigh16 = vmulq_n_u16(Ghigh16,(uint16_t)(gKernel_far[iter]*64));
            Rhigh16 = vmulq_n_u16(Rhigh16,(uint16_t)(gKernel_far[iter]*64));


            Bhigh16 = vshrq_n_u16(Bhigh16,6);
            Ghigh16 = vshrq_n_u16(Ghigh16,6);
            Rhigh16 = vshrq_n_u16(Rhigh16,6);

            Bhigh = vqmovn_u16(Bhigh16);
            Ghigh = vqmovn_u16(Ghigh16);
            Rhigh = vqmovn_u16(Rhigh16);

            Bvector = vcombine_u8(Blow,Bhigh);
            Gvector = vcombine_u8(Glow,Ghigh);
            Rvector = vcombine_u8(Rlow,Rhigh);



            pixelChannels.val[0] = Bvector;
            pixelChannels.val[1] = Gvector;
            pixelChannels.val[2] = Rvector;
            pixelChannels.val[3] = Avector;


            vst4q_u8(arrayQptr,pixelChannels);

            arrayInPtr += 64;
//            arrayOutPtr += 16;
            arrayQptr += 64;
        }

        for(int j=0;j<=height;j++)
        {
            for(int i=0;i<=width;i++)
            {
                if((i-ceil_rFar+iter)>=0 && (i-ceil_rFar+iter) < width)
                {
                    arrayOutPtr[(j*width+i)*4+0] += q[(j*width+i-ceil_rFar+iter)*4+0];
                    arrayOutPtr[(j*width+i)*4+1] += q[(j*width+i-ceil_rFar+iter)*4+1];
                    arrayOutPtr[(j*width+i)*4+2] += q[(j*width+i-ceil_rFar+iter)*4+2];
                    arrayOutPtr[(j*width+i)*4+3] = q[(j*width+i-ceil_rFar+iter)*4+3];
                }
            }
        }

    }


    for(int iter=0; iter<=2*ceil_rFar; iter++)
    {
        for (int i=0; i <width*height; i+= 16){
            uint8x16x4_t pixelChannels = vld4q_u8(arrayInPtr);
            int y = i/width;

            uint8x16_t Bvector = pixelChannels.val[0];
            uint8x16_t Gvector = pixelChannels.val[1];
            uint8x16_t Rvector = pixelChannels.val[2];
            uint8x16_t Avector = pixelChannels.val[3];

            uint8x8_t Blow = vget_low_u8(Bvector);
            uint8x8_t Glow = vget_low_u8(Gvector);
            uint8x8_t Rlow = vget_low_u8(Rvector);

            uint8x8_t Bhigh = vget_high_u8(Bvector);
            uint8x8_t Ghigh = vget_high_u8(Gvector);
            uint8x8_t Rhigh = vget_high_u8(Rvector);

            uint16x8_t Bhigh16 = vmovl_u8(Bhigh);
            uint16x8_t Ghigh16 = vmovl_u8(Ghigh);
            uint16x8_t Rhigh16 = vmovl_u8(Rhigh);

            uint16x8_t Blow16 = vmovl_u8(Blow);
            uint16x8_t Glow16 = vmovl_u8(Glow);
            uint16x8_t Rlow16 = vmovl_u8(Rlow);

            if(y<a0)
            {
                if(sigma_far > 0.6)
                {
                    Blow16 = vmulq_n_u16(Blow16,(uint16_t)(gKernel_far[iter]*64));
                    Glow16 = vmulq_n_u16(Glow16,(uint16_t)(gKernel_far[iter]*64));
                    Rlow16 = vmulq_n_u16(Rlow16,(uint16_t)(gKernel_far[iter]*64));
                }
            }
            else if(y >= a0 && y < a1)
            {
                float sigmaFar_2 = (a1-y)/(a1-a0);
                sigmaFar_2 = sigma_far*sigmaFar_2;
                if(sigmaFar_2 > 0.6)
                {
                    gKernel_create_vector(sigmaFar_2,2*ceil_rFar,gKernel_far_two);
                }
                Blow16 = vmulq_n_u16(Blow16,(uint16_t)(gKernel_far[iter]*64));
                Glow16 = vmulq_n_u16(Glow16,(uint16_t)(gKernel_far[iter]*64));
                Rlow16 = vmulq_n_u16(Rlow16,(uint16_t)(gKernel_far[iter]*64));
            }
            else if(y>=a1 && y < a2)
            {
                goto here;
            }
            else if(y >= a2 && y < a3)
            {
                float sigmaNear_2 = (y-a2)/(a2-a1);
                sigmaNear_2 = sigmaNear_2*sigma_near;
                if(sigmaNear_2 > 0.6)
                {
                    gKernel_create_vector(sigmaNear_2,2*ceil_rNear,gKernel_near_two);
                }
                Blow16 = vmulq_n_u16(Blow16,(uint16_t)(gKernel_near[iter]*64));
                Glow16 = vmulq_n_u16(Glow16,(uint16_t)(gKernel_near[iter]*64));
                Rlow16 = vmulq_n_u16(Rlow16,(uint16_t)(gKernel_near[iter]*64));
            }
            else
            {
                if(sigma_near > 0.6)
                {
                    Blow16 = vmulq_n_u16(Blow16,(uint16_t)(gKernel_near[iter]*64));
                    Glow16 = vmulq_n_u16(Glow16,(uint16_t)(gKernel_near[iter]*64));
                    Rlow16 = vmulq_n_u16(Rlow16,(uint16_t)(gKernel_near[iter]*64));
                }

            }




            Blow16 = vshrq_n_u16(Blow16,6);
            Glow16 = vshrq_n_u16(Rlow16,6);
            Rlow16 = vshrq_n_u16(Rlow16,6);

            Blow = vqmovn_u16(Blow16);
            Glow = vqmovn_u16(Glow16);
            Rlow = vqmovn_u16(Rlow16);



            if(y<a0)
            {
                if(sigma_far > 0.6)
                {
                    Bhigh16 = vmulq_n_u16(Bhigh16,(uint16_t)(gKernel_far[iter]*64));
                    Ghigh16 = vmulq_n_u16(Ghigh16,(uint16_t)(gKernel_far[iter]*64));
                    Rhigh16 = vmulq_n_u16(Rhigh16,(uint16_t)(gKernel_far[iter]*64));
                }
            }
            else if(y >= a0 && y < a1)
            {
                float sigmaFar_2 = (a1-y)/(a1-a0);
                sigmaFar_2 = sigma_far*sigmaFar_2;
                if(sigmaFar_2 > 0.6)
                {
                    gKernel_create_vector(sigmaFar_2,2*ceil_rFar,gKernel_far_two);
                }
                Bhigh16 = vmulq_n_u16(Bhigh16,(uint16_t)(gKernel_far[iter]*64));
                Ghigh16 = vmulq_n_u16(Ghigh16,(uint16_t)(gKernel_far[iter]*64));
                Rhigh16 = vmulq_n_u16(Rhigh16,(uint16_t)(gKernel_far[iter]*64));
            }
            else if(y >= a2 && y < a3)
            {
                float sigmaNear_2 = (y-a2)/(a2-a1);
                sigmaNear_2 = sigmaNear_2*sigma_near;
                if(sigmaNear_2 > 0.6)
                {
                    gKernel_create_vector(sigmaNear_2,2*ceil_rNear,gKernel_near_two);
                }
                Bhigh16 = vmulq_n_u16(Bhigh16,(uint16_t)(gKernel_near[iter]*64));
                Ghigh16 = vmulq_n_u16(Ghigh16,(uint16_t)(gKernel_near[iter]*64));
                Rhigh16 = vmulq_n_u16(Rhigh16,(uint16_t)(gKernel_near[iter]*64));
            }
            else
            {
                if(sigma_near > 0.6)
                {
                    Bhigh16 = vmulq_n_u16(Bhigh16,(uint16_t)(gKernel_near[iter]*64));
                    Ghigh16 = vmulq_n_u16(Ghigh16,(uint16_t)(gKernel_near[iter]*64));
                    Rhigh16 = vmulq_n_u16(Rhigh16,(uint16_t)(gKernel_near[iter]*64));
                }

            }



            Bhigh16 = vshrq_n_u16(Bhigh16,6);
            Ghigh16 = vshrq_n_u16(Ghigh16,6);
            Rhigh16 = vshrq_n_u16(Rhigh16,6);

            Bhigh = vqmovn_u16(Bhigh16);
            Ghigh = vqmovn_u16(Ghigh16);
            Rhigh = vqmovn_u16(Rhigh16);

            Bvector = vcombine_u8(Blow,Bhigh);
            Gvector = vcombine_u8(Glow,Ghigh);
            Rvector = vcombine_u8(Rlow,Rhigh);



            pixelChannels.val[0] = Bvector;
            pixelChannels.val[1] = Gvector;
            pixelChannels.val[2] = Rvector;
            pixelChannels.val[3] = Avector;

            here:
            vst4q_u8(arrayQptr,pixelChannels);

            arrayInPtr += 64;
//            arrayOutPtr += 16;
            arrayQptr += 64;
        }

        for(int j=0;j<=height;j++)
        {
            for(int i=0;i<=width;i++)
            {
                if((i-ceil_rFar+iter)>=0 && (i-ceil_rFar+iter) < width)
                {
                    arrayOutPtr[(j*width+i)*4+0] += q[(j-ceil_rFar+iter*width+i)*4+0];
                    arrayOutPtr[(j*width+i)*4+1] += q[(j-ceil_rFar+iter*width+i)*4+1];
                    arrayOutPtr[(j*width+i)*4+2] += q[(j-ceil_rFar+iter*width+i)*4+2];
                    arrayOutPtr[(j*width+i)*4+3] = q[(j-ceil_rFar+iter*width+i)*4+3];
                }
            }
        }

    }



//    for (int j=0;j<height;j++){
//        for (int i=0;i<width;i++) {
//            int B = pixels[j*width+i]%0x100;
//            int G = (pixels[j*width+i]>>8)%0x100;
//            int R = (pixels[j*width+i]>>16)%0x100;
//            int A = 0xff;
//            R=0;
//            int color = (A & 0xff) << 24 | (R & 0xff) << 16 | (G & 0xff) << 8 | (B & 0xff);
//
//            outputPixels[j*width+i]=color;
//        }
//    }

    env->ReleaseIntArrayElements(inputPixels_, pixels, 0);
    env->ReleaseIntArrayElements(outputPixels_, outputPixels, 0);
    return 0;
}