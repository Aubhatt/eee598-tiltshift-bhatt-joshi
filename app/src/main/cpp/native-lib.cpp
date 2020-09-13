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

// TODO: Parallelize using std::for_each and vectors

jfloat* gaussian1D_kernel(jfloat sigma, jint k_radius) {
    int k_width = 2*k_radius + 1;

    int x;
    jfloat r, s = 2.0 * sigma * sigma;
    jfloat sum=0;

    // Filling kernel elements
    auto *kernel = new jfloat[k_width];

    for(int i=0; i<k_width; i++) {
        x = i-k_radius;
        kernel[i] = (exp(-(x*x) / s)) / (M_PI * s);
        sum += kernel[i];
    }

    // Normalize weights (so that the sums add to 1)
    for(int i=0; i<k_width; i++) {
        kernel[i] /= sum;
    }

    return kernel;
}

jfloat* gaussian_kernel(jfloat sigma, jint k_radius) {
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
                  jint height,
                  jint k_radius) {
    x_start = fmax(x_start, 0);
    y_start = fmax(y_start, 0);

    x_end = fmin(x_end, width);
    y_end = fmin(y_end, width);

    // Convolution of image with kernel a.k.a. applying filter
    for (int j=y_start; j<y_end; j++){
        for (int i=x_start + k_radius; i<x_end - k_radius; i++) {
            outputPixels[j*width+i] = pixels[j*width + i];
        }
    }
}

void apply_filterFast(jint *pixels,
                  jfloat *kernel,
                  jint *outputPixels,
                  jint x_start,
                  jint y_start,
                  jint x_end,
                  jint y_end,
                  jint width,
                  jint height,
                  jint k_radius) {

    auto* tempPixels = new jint[width*height];
    int k_width = 2*k_radius + 1;

    x_start = fmax(x_start, k_radius);
    y_start = fmax(y_start, k_radius);

    x_end = fmin(x_end, width-k_radius);
    y_end = fmin(y_end, height-k_radius);

    // Convolution of image with kernel a.k.a. applying filter

    // First pass
    for (int j=y_start; j<y_end; j++){
        for (int i=x_start; i<x_end; i++) {
            jfloat B = 0, G = 0, R = 0, A = 0xff;

            //Applying kernel vertically to pixel
            for(int k_y=0; k_y<k_width; k_y++) {
                int y = k_y + j - k_radius;
                int x = i;

                uint32_t b = pixels[y*width + x] & 0xFF; //% 0x100;
                uint32_t g = (pixels[y*width + x] >> 8) & 0xFF;
                uint32_t r = (pixels[y*width + x] >> 16) & 0xFF;

                B +=  (b*kernel[k_y]);
                G +=  (g*kernel[k_y]);
                R +=  (r*kernel[k_y]);
            }

            uint32_t _B = (uint32_t) B;
            uint32_t _G = (uint32_t) G;
            uint32_t _R = (uint32_t) R;
            uint32_t _A = (uint32_t) A;
            uint32_t color = (_A & 0xff) << 24 | (_R & 0xff) << 16 | (_G & 0xff) << 8 | (_B & 0xff);

            tempPixels[j*width+i]=color;
        }
    }

    // Second pass
    for (int j=y_start; j<y_end; j++){
        for (int i=x_start; i<x_end; i++) {
            jfloat B = 0, G = 0, R = 0, A = 0xff;

            //Applying kernel horizontally to pixel
            for(int k_x=0; k_x<k_width; k_x++) {
                int x = k_x + i - k_radius;
                int y = j;

                uint32_t b = tempPixels[y*width + x] & 0xFF; //% 0x100;
                uint32_t g = (tempPixels[y*width + x] >> 8) & 0xFF;
                uint32_t r = (tempPixels[y*width + x] >> 16) & 0xFF;

                B +=  (b*kernel[k_x]);
                G +=  (g*kernel[k_x]);
                R +=  (r*kernel[k_x]);
            }

            uint32_t _B = (uint32_t) B;
            uint32_t _G = (uint32_t) G;
            uint32_t _R = (uint32_t) R;
            uint32_t _A = (uint32_t) A;
            uint32_t color = (_A & 0xff) << 24 | (_R & 0xff) << 16 | (_G & 0xff) << 8 | (_B & 0xff);

            outputPixels[j*width+i]=color;
        }
    }

    delete[] tempPixels; // Free intermediate pixel storage
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
    y_end = fmin(y_end, height-k_radius);

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

void gaussian_filter(jint *pixels,
                             jint *outputPixels,
                             jint x_start,
                             jint y_start,
                             jint x_end,
                             jint y_end,
                             jint width,
                             jint height,
                             jint sigma,
                             jint k_radius,
                             jint fast) {

    jfloat *kernel;

    // Only apply if sigma is greater than 0.6
    if(sigma >= 0.6) {
        if(fast) {
            kernel = gaussian1D_kernel(sigma, k_radius);
            apply_filterFast(pixels, kernel, outputPixels, x_start, y_start, x_end, y_end, width, height, k_radius);
        }
        else {
            kernel = gaussian_kernel(sigma, k_radius);
            apply_filter(pixels, kernel, outputPixels, x_start, y_start, x_end, y_end, width, height, k_radius);
        }
        delete[] kernel; // Free allocated kernel
    }
    else {
        copy_buffer2D(pixels, outputPixels, x_start, y_start, x_end, y_end, width, height, k_radius);
    }

}

void gaussianGradient_filter(jint *pixels,
        jint *outputPixels,
        jint x_start,
        jint y_start,
        jint x_end,
        jint y_end,
        jint width,
        jint height,
        jint sigma,
        jint k_radius,
        jint climb,
        jint fast) { // climb: 1 - no_blur to blur; 0 - blur to no_blur

    jfloat *kernel;
    jfloat sigma_grad;

    // Apply different kernel for every row in the image to generate a gradient blur
    for(int j=y_start; j<y_end; j++) {

        // Check if the blur is increasing or decreasing
        if(climb)
            sigma_grad = sigma * abs(j - y_start) / abs(y_end - y_start);
        else
            sigma_grad = sigma * abs(j - y_end) / abs(y_end - y_start);

        // Only apply if sigma is greater than 0.6
        if(sigma_grad >= 0.6) {
            if(fast) {
                kernel = gaussian1D_kernel(sigma_grad, k_radius);
                apply_filterFast(pixels, kernel, outputPixels, x_start, j, x_end, j+1, width, height, k_radius);
            }
            else {
                kernel = gaussian_kernel(sigma_grad, k_radius);
                apply_filter(pixels, kernel, outputPixels, x_start, j, x_end, j+1, width, height, k_radius);
            }
            delete[] kernel; // Free allocated kernel
        }
        else {
            copy_buffer2D(pixels, outputPixels, x_start, j, x_end, j+1, width, height, k_radius);
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

    jint *pixels = env->GetIntArrayElements(inputPixels_, NULL);
    jint *outputPixels = env->GetIntArrayElements(outputPixels_, NULL);

    // Selecting convolution algorithm
    jint fast = 1; // 1: Weight Vector, 0: Weight Matrix

    // Selecting kernel radius as per the sigma values
    jint k_radius = ceil(2*fmax(sigma_far, sigma_near));

    // Create threads for each strip
    std::thread strip0 (gaussian_filter, pixels, outputPixels, 0,0, width, a0, width, height, sigma_far, k_radius, fast);
    std::thread strip1 (gaussianGradient_filter, pixels, outputPixels, 0, a0, width, a1, width, height, sigma_far, k_radius, 0, fast);
    std::thread strip2 (gaussian_filter, pixels, outputPixels, 0, a1, width, a2, width, height, 0, k_radius, fast);
    std::thread strip3 (gaussianGradient_filter, pixels, outputPixels, 0, a2, width, a3, width, height, sigma_far, k_radius, 1, fast);
    std::thread strip4 (gaussian_filter, pixels, outputPixels, 0, a3, width, height, width, height, sigma_far, k_radius, fast);

    // Wait for the threads to finish
    strip0.join();
    strip1.join();
    strip2.join();
    strip3.join();
    strip4.join();

    // Release the pixels
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