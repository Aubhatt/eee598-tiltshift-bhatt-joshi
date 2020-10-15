package edu.asu.ame.meteor.speedytiltshift2018;

import android.util.Log;

public class SpeedyTiltShiftJava extends Thread{
    private int[] pixels;
    private int[] pixelsOut;
    private int x_start;
    private int y_start;
    private int x_end;
    private int y_end;
    private int width;
    private int height;
    private float sigma;
    private int k_radius;
    private int climb;
    private int fast;

    private static double[] gaussian_kernel(float sigma, int k_radius) {
        int k_height = 2*k_radius + 1;
        int k_width = 2*k_radius + 1;

        int x, y;
        double r, s = (2.0 * sigma * sigma);

        double sum= 0.0;

        // Filling kernel elements
        double[] kernel = new double[k_height*k_width];
        for(int j=0; j<k_height; j++) {
            for(int i=0; i<k_width; i++) {
                x = i-k_radius;
                y = j-k_radius;
                r = Math.sqrt(x*x + y*y);
                kernel[j*k_width + i] = ((Math.exp(-(r*r) / s)) / (Math.PI * s));
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

    private static void copy_buffer2D(final int[] pixels,
                                      int[] pixelsOut,
                                      int x_start,
                                      int y_start,
                                      int x_end,
                                      int y_end,
                                      int width,
                                      int height,
                                      int k_radius) {
        x_start = Integer.max(x_start, 0);
        y_start = Integer.max(y_start, 0);

        x_end = Integer.max(x_end, width);
        y_end = Integer.max(y_end, height);

        // Convolution of image with kernel a.k.a. applying filter
        for (int j=y_start; j<y_end; j++){
            for (int i=x_start; i<x_end; i++) {
                pixelsOut[j*width+i] = pixels[j*width + i];
            }
        }
    }

    private static void apply_filter(final int[] pixels,
                                     final double[] kernel,
                                     int[] pixelsOut,
                                     int x_start,
                                     int y_start,
                                     int x_end,
                                     int y_end,
                                     int width,
                                     int height,
                                     int k_radius) {
        double B, G, R;
        int x, y;
        int b, g, r;
        int _B, _G, _R, _A;
        int color;

        int k_height = 2*k_radius + 1;
        int k_width = 2*k_radius + 1;

        x_start = Integer.max(x_start, 0);
        y_start = Integer.max(y_start, 0);

        x_end = Integer.min(x_end, width);
        y_end = Integer.min(y_end, height);

        // Convolution of image with kernel a.k.a. applying filter
        for (int j=y_start; j<y_end; j++){
            for(int i=x_start; i<x_end; i++) {
                B = 0; G = 0; R = 0; _A = 0xff;

                //Apply kernel to pixel
                for(int k_y=0; k_y<k_height; k_y++) {
                    for(int k_x=0; k_x<k_width; k_x++) {
                        y = k_y + j - k_radius;
                        x = k_x + i - k_radius;

                        if( (y >= 0) && (y < height) && (x >= 0) && (x < width) ) {
                            b = pixels[y*width + x] & 0xFF; //% 0x100;
                            g = (pixels[y*width + x] >> 8) & 0xFF;
                            r = (pixels[y*width + x] >> 16) & 0xFF;

                            //kernel[k_y*k_width + k_x] = 1;
                            B +=  b*kernel[k_y*k_width + k_x];
                            G +=  g*kernel[k_y*k_width + k_x];
                            R +=  r*kernel[k_y*k_width + k_x];
                        }
                    }
                }

                _B = (int) B;
                _G = (int) G;
                _R = (int) R;
                color = (_A & 0xff) << 24 | (_R & 0xff) << 16 | (_G & 0xff) << 8 | (_B & 0xff);
                pixelsOut[j*width + i]=color;
            }
        }
    }

    private static void gaussian_filter(int[] pixels,
                                        int[] pixelsOut,
                                        int x_start,
                                        int y_start,
                                        int x_end,
                                        int y_end,
                                        int width,
                                        int height,
                                        float sigma,
                                        int k_radius,
                                        int fast) {
        double[] kernel;

        // Only apply if sigma is greater than 0.6
        if(sigma >= 0.6) {
            if(fast == 1) {
                kernel = gaussian_kernel(sigma, k_radius);
                apply_filter(pixels, kernel, pixelsOut, x_start, y_start, x_end, y_end, width, height, k_radius);
            }
            else {
                kernel = gaussian_kernel(sigma, k_radius);
                apply_filter(pixels, kernel, pixelsOut, x_start, y_start, x_end, y_end, width, height, k_radius);
            }
        }
        else {
            copy_buffer2D(pixels, pixelsOut, x_start, y_start, x_end, y_end, width, height, k_radius);
        }
    }

    private static void gaussianGradient_filter(int[] pixels,
                                                int[] pixelsOut,
                                                int x_start,
                                                int y_start,
                                                int x_end,
                                                int y_end,
                                                int width,
                                                int height,
                                                float sigma,
                                                int k_radius,
                                                int climb,
                                                int fast) {
        double[] kernel;
        float sigma_grad;

        for(int j=y_start; j<y_end; j++) {
            if(climb == 1)
                sigma_grad = (float)((double)sigma * (float)Math.abs(j-y_start) / (float)Math.abs(y_end - y_start));
            else
                sigma_grad = (float)((double)sigma * (float)Math.abs(j-y_end) / (float)Math.abs(y_end - y_start));

            // Only apply if sigma is greater than 0.6
            if(sigma_grad >= 0.6) {
                if(fast == 1) {
                    kernel = gaussian_kernel(sigma_grad, k_radius);
                    apply_filter(pixels, kernel, pixelsOut, x_start, j, x_end, j+1, width, height, k_radius);
                }
                else {
                    kernel = gaussian_kernel(sigma_grad, k_radius);
                    apply_filter(pixels, kernel, pixelsOut, x_start, j, x_end, j+1, width, height, k_radius);
                }
            }
            else {
                copy_buffer2D(pixels, pixelsOut, x_start, j, x_end, j+1, width, height, k_radius);
            }
        }

    }

    SpeedyTiltShiftJava(int[] pixels,
                        int[] pixelsOut,
                        int x_start,
                        int y_start,
                        int x_end,
                        int y_end,
                        int width,
                        int height,
                        float sigma,
                        int k_radius,
                        int climb,
                        int fast) {

        this.pixels =  pixels;
        this.pixelsOut = pixelsOut;
        this.x_start = x_start;
        this.y_start = y_start;
        this.x_end = x_end;
        this.y_end = y_end;
        this.width = width;
        this.height = height;
        this.sigma = sigma;
        this.k_radius = k_radius;
        this.climb = climb;
        this.fast = fast;
    }

    public void run() {

        if(this.climb == 2) {
            gaussian_filter(this.pixels,
                    this.pixelsOut,
                    this.x_start,
                    this.y_start,
                    this.x_end,
                    this.y_end,
                    this.width,
                    this.height,
                    this.sigma,
                    this.k_radius,
                    this.fast);
        }else {
            gaussianGradient_filter(this.pixels,
                    this.pixelsOut,
                    this.x_start,
                    this.y_start,
                    this.x_end,
                    this.y_end,
                    this.width,
                    this.height,
                    this.sigma,
                    this.k_radius,
                    this.climb,
                    this.fast);
        }
    }

    public void printInfo(int a) {
        Log.d("SPEEDY_TS", "printInfo" + String.valueOf(a) +String.valueOf(this.x_start) + ", " + String.valueOf(this.y_start) + ", " + String.valueOf(this.x_end) + ", " + String.valueOf(this.y_end));
    }
}