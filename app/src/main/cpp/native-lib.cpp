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

jfloat* gaussianGradient_filter(jfloat sigma, jint k_radius, jint y1, jint y2) {
    int k_height = 2*k_radius + 1;
    int k_width = 2*k_radius + 1;

    int x, y;
    jfloat r, s, sigma_grad;

    jfloat sum=0;

    // Filling kernel elements
    auto *kernel = new jfloat[k_height*k_width];
    for(int j=0; j<k_height; j++) {
        for(int i=0; i<k_width; i++) {
            x = i-k_radius;
            y = j-k_radius;
            sigma_grad = sigma * abs(y1 - j) / abs(y1 - y2);
            r = sqrt(x*x + y*y);
            s = 2.0 * sigma_grad * sigma_grad;
            kernel[j*k_width + i] = (exp(-(r*r) / s)) / (M_PI * s);
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

jfloat* gaussian_filter(jfloat sigma, jint k_radius) {
    int k_height = 2*k_radius + 1;
    int k_width = 2*k_radius + 1;

    int x, y;
    jfloat r, s = 2.0 * sigma * sigma;

    jfloat sum=0;

    // Filling kernel elements
    auto *kernel = new jfloat[k_height*k_width];
    for(int j=0; j<k_height; j++) {
        for(int i=0; i<k_width; i++) {
            x = i-k_radius;
            y = j-k_radius;
            r = sqrt(x*x + y*y);
            kernel[j*k_width + i] = (exp(-(r*r) / s)) / (M_PI * s);
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

void copy_buffer2D(jint *pixels,
                  jint *outputPixels,
                  jint x_start,
                  jint y_start,
                  jint x_end,
                  jint y_end,
                  jint width,
                  jint height) {
    x_start = fmax(x_start, 0);
    y_start = fmax(y_start, 0);

    x_end = fmin(x_end, width);
    y_end = fmin(y_end, width);

    // Convolution of image with kernel a.k.a. applying filter
    for (int j=y_start; j<y_end; j++){
        for (int i=x_start; i<x_end; i++) {
            outputPixels[j*width+i] = pixels[j*width + i];
        }
    }
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
            jfloat B = 0, G = 0, R = 0, A = 0xff;

            //Applying kernel to pixel
            for(int k_y=0; k_y<k_height; k_y++) {
                for(int k_x=0; k_x<k_width; k_x++) {
                    int y = k_y + j - k_radius;
                    int x = k_x + i - k_radius;

                    uint32_t b = pixels[y*width + x] & 0xFF; //% 0x100;
                    uint32_t g = (pixels[y*width + x] >> 8) & 0xFF;
                    uint32_t r = (pixels[y*width + x] >> 16) & 0xFF;

                    B +=  (b*kernel[k_y*k_width + k_x]);
                    G +=  (g*kernel[k_y*k_width + k_x]);
                    R +=  (r*kernel[k_y*k_width + k_x]);
                }
            }
            uint32_t _B = (uint32_t) B;
            uint32_t _G = (uint32_t) G;
            uint32_t _R = (uint32_t) R;
            uint32_t _A = (uint32_t) A;
            uint32_t color = (_A & 0xff) << 24 | (_R & 0xff) << 16 | (_G & 0xff) << 8 | (_B & 0xff);

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

    LOGI("Sigma_far = %f", sigma_far);

    // Creating filters
    jint k_radius0 = ceil(2*sigma_far);
    jint k_radius1 = ceil(2*sigma_far/2);
    jint k_radius2 = ceil(2*sigma_near/2);
    jint k_radius3 = ceil(2*sigma_near);

    jfloat *kernel_l0 = gaussian_filter(sigma_far, k_radius0);
    jfloat *kernel_l1 = gaussianGradient_filter(sigma_far, k_radius1, a1, a0);
    jfloat *kernel_l2 = gaussianGradient_filter(sigma_near, k_radius2, a2, a3);
    jfloat *kernel_l3 = gaussian_filter(sigma_near, k_radius3);

    auto *kernel_dummy = new jfloat[1];
    kernel_dummy[0] = 1;

    // Applying the filters on the input image
    apply_filter(pixels, kernel_l0, outputPixels,0,0, width, a0, width, height, k_radius0);
    apply_filter(pixels, kernel_l1, outputPixels,0,a0, width, a1, width, height, k_radius1);
    copy_buffer2D(pixels, outputPixels, 0, a1, width, a2, width, height);
    apply_filter(pixels, kernel_l2, outputPixels,0,a2, width, a3, width, height, k_radius2);
    apply_filter(pixels, kernel_l3, outputPixels,0,a3, width, height, width, height, k_radius3);

//     std::thread filter (apply_filter, pixels, kernel_l0, outputPixels, 0, 0, width, height, width, height, k_radius0);
//     filter.join();

     LOGI("Done");

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