#include <jni.h>
#include <string>
#include <cpu-features.h>
#include <android/log.h>
#include <thread>
#include <cmath>
#include <algorithm>

#define TAG "SPEEDY_TS"

#define LOGE(...) __android_log_print(ANDROID_LOG_ERROR, TAG, __VA_ARGS__)
#define LOGW(...) __android_log_print(ANDROID_LOG_WARN, TAG, __VA_ARGS__)
#define LOGI(...) __android_log_print(ANDROID_LOG_INFO, TAG, __VA_ARGS__)
#define LOGD(...) __android_log_print(ANDROID_LOG_DEBUG, TAG, __VA_ARGS__)

jfloat* avg_filter(jint k_radius) {
    int k_height = 2*k_radius + 1;
    int k_width = 2*k_radius + 1;

    // Filling kernel elements
    auto *kernel = new jfloat[k_height*k_width];
    for(int j=0; j<k_height; j++) {
        for(int i=0; i<k_width; i++) {
            kernel[j*k_width + i] = 1.00/(k_height*k_width);
        }
    }
    return kernel;
}

jfloat* gaussian_filter(jint k_radius,
        jfloat sigma) {
    int k_height = 2*k_radius + 1;
    int k_width = 2*k_radius + 1;

    jfloat sum=0;

    // Filling kernel elements
    auto *kernel = new jfloat[k_height*k_width];
    for(int j=0; j<k_height; j++) {
        for(int i=0; i<k_width; i++) {
            kernel[j*k_width + i] = pow(M_E, (-(pow(i-k_radius, 2) + pow(j-k_radius, 2))/(2*pow(sigma,2)))) / (2*M_PI*pow(sigma, 2));
            sum += kernel[j*k_width + i];
        }
    }

    // Normalize weights (so that the sums add to 1)
    for(int j=0; j<k_height; j++) {
        for(int i=0; i<k_width; i++) {
            kernel[j*k_width + i] /= sum;
        }
    }

    return kernel;
}

void apply_filter(jint *pixels,
        jfloat *kernel,
        jint *outputPixels,
        jint x_start,
        jint y_start,
        jint x_end,
        jint y_end,
        jint width,
        jint height,
        jint k_radius) {

    int k_height = 2*k_radius + 1;
    int k_width = 2*k_radius + 1;

    x_start = fmax(x_start, k_radius);
    y_start = fmax(y_start, k_radius);

    x_end = fmin(x_end, width-k_radius);
    y_end = fmin(y_end, width-k_radius);

    // Convolution of image with kernel a.k.a. applying filter
    for (int j=y_start; j<y_end; j++){
        for (int i=x_start; i<x_end; i++) {
            uint8_t B = 0, G = 0, R = 0, A = 0xff;

            //Applying kernel to pixel
            for(int k_y=0; k_y<k_height; k_y++) {
                for(int k_x=0; k_x<k_width; k_x++) {
                    int y = k_y + j - k_radius;
                    int x = k_x + i - k_radius;

                    uint8_t b = pixels[y*width + x] % 0x100;
                    uint8_t g = (pixels[y*width + x] >> 8) % 0x100;
                    uint8_t r = (pixels[y*width + x] >> 16) % 0x100;

                    B += (uint8_t) (b*kernel[k_y*k_width + k_x]);
                    G += (uint8_t) (g*kernel[k_y*k_width + k_x]);
                    R += (uint8_t) (r*kernel[k_y*k_width + k_x]);
                }
            }

            uint32_t color = (A & 0xff) << 24 | (R & 0xff) << 16 | (G & 0xff) << 8 | (B & 0xff);

            outputPixels[j*width+i]=color;
        }
    }
}

extern "C"
JNIEXPORT jint JNICALL
Java_edu_asu_ame_meteor_speedytiltshift2018_SpeedyTiltShift_tiltshiftcppnative(JNIEnv *env,
                                                                               jobject instance,
                                                                               jintArray inputPixels_,
                                                                               jintArray outputPixels_,
                                                                               jint width,
                                                                               jint height,
                                                                               jfloat sigma_far,
                                                                               jfloat sigma_near,
                                                                               jint a0, jint a1,
                                                                               jint a2, jint a3) {

    /*
     * a3 a2 a1 a0 [JAVA] [C++] [NEON] sigma_far sigma_near
     */

    jint *pixels = env->GetIntArrayElements(inputPixels_, NULL);
    jint *outputPixels = env->GetIntArrayElements(outputPixels_, NULL);

    // Create filter
    jint k_radius = 2;
    jfloat *kernel = gaussian_filter(k_radius, 2);

    // Applying gaussian filter on the input image
    apply_filter(pixels, kernel, outputPixels, 0, 0, width/2, height/2, width, height, k_radius);
    apply_filter(pixels, kernel, outputPixels, width/2, height/2, width, height, width, height, k_radius);
    // std::thread filter (apply_filter, pixels, kernel, outputPixels, 0, 0, width/2, height/2, width, height, k_radius);
    // filter.join();

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

    for (int j=0;j<height;j++){
        for (int i=0;i<width;i++) {
            int B = pixels[j*width+i]%0x100;
            int G = (pixels[j*width+i]>>8)%0x100;
            int R = (pixels[j*width+i]>>16)%0x100;
            int A = 0xff;
            R=0;
            int color = (A & 0xff) << 24 | (R & 0xff) << 16 | (G & 0xff) << 8 | (B & 0xff);

            outputPixels[j*width+i]=color;
        }
    }

    env->ReleaseIntArrayElements(inputPixels_, pixels, 0);
    env->ReleaseIntArrayElements(outputPixels_, outputPixels, 0);
    return 0;
}