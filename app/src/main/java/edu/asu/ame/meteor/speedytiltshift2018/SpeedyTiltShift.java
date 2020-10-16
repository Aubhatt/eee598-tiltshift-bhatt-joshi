package edu.asu.ame.meteor.speedytiltshift2018;

import android.graphics.Bitmap;
import android.util.Log;

import java.util.Arrays;

import static java.lang.Float.max;
import static java.lang.Math.ceil;

public class SpeedyTiltShift {
    static SpeedyTiltShift Singleton = new SpeedyTiltShift();

    // Used to load the 'native-lib' library on application startup.
    static {
        System.loadLibrary("native-lib");
    }

    private static double[] gaussian1D_kernel(float sigma, int k_radius) {
        int k_width = 2*k_radius + 1;
        double r, s = (2.0 * sigma * sigma);
        double sum= 0.0;
        double x;

        // Filling kernel elements
        double[] kernel = new double[k_width];

        for(int i=0; i<k_width; i++) {
            x = i-k_radius;
            kernel[i] = ((Math.exp(-(x*x) / s)) / (Math.PI * s));
            sum += kernel[i];
        }

        // Normalize weights (so that the sums add to 1)
        for(int i=0; i<k_width; i++) {
            kernel[i] /= sum;
        }

        return kernel;
    }


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

    private static void apply_filterFast(final int[] pixels,
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

        int [] tempPixels = new int[width*height];
        int k_width = 2*k_radius + 1;

        x_start = Integer.max(x_start, 0);
        y_start = Integer.max(y_start, 0);

        x_end = Integer.min(x_end, width);
        y_end = Integer.min(y_end, height);

        // First pass
        for (int j=y_start; j<y_end; j++){
            for (int i=x_start; i<x_end; i++) {
                B = 0; G = 0; R = 0; _A = 0xff;

                //Applying kernel vertically to pixel
                for(int k_y=0; k_y<k_width; k_y++) {
                    y = k_y + j - k_radius;
                    x = i;

                    if( (y >= 0) && (y < height) ) {
                        b = pixels[y * width + x] & 0xFF; //% 0x100;
                        g = (pixels[y * width + x] >> 8) & 0xFF;
                        r = (pixels[y * width + x] >> 16) & 0xFF;

                        // Notice that B, G and R are float and not int
                        // to avoid accumulating rounding error
                        B += (b * kernel[k_y]);
                        G += (g * kernel[k_y]);
                        R += (r * kernel[k_y]);
                    }
                }

                // Now we convert the accumulated float values to int values
                _B = (int) B;
                _G = (int) G;
                _R = (int) R;
                color = (_A & 0xff) << 24 | (_R & 0xff) << 16 | (_G & 0xff) << 8 | (_B & 0xff);

                tempPixels[j*width+i]=color;
            }
        }

        for (int j=y_start; j<y_end; j++){
            for (int i=x_start; i<x_end; i++) {
                B = 0; G = 0; R = 0; _A = 0xff;

                //Applying kernel horizontally to pixel
                for(int k_x=0; k_x<k_width; k_x++) {
                    x = k_x + i - k_radius;
                    y = j;

                    if( (x >= 0) && (x < width) ) {
                        b = tempPixels[y * width + x] & 0xFF; //% 0x100;
                        g = (tempPixels[y * width + x] >> 8) & 0xFF;
                        r = (tempPixels[y * width + x] >> 16) & 0xFF;



                        // Notice that B, G and R are float and not int
                        // to avoid accumulating rounding error
                        B += (b * kernel[k_x]);
                        G += (g * kernel[k_x]);
                        R += (r * kernel[k_x]);
                    }
                }

                // Now we convert the accumulated float values to int values
                _B = (int) B;
                _G = (int) G;
                _R = (int) R;
                color = (_A & 0xff) << 24 | (_R & 0xff) << 16 | (_G & 0xff) << 8 | (_B & 0xff);

                pixelsOut[j*width+i]=color;
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
                kernel = gaussian1D_kernel(sigma, k_radius);
                apply_filterFast(pixels, kernel, pixelsOut, x_start, y_start, x_end, y_end, width, height, k_radius);
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
                    kernel = gaussian1D_kernel(sigma_grad, k_radius);
                    apply_filterFast(pixels, kernel, pixelsOut, x_start, j, x_end, j+1, width, height, k_radius);
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

    public static Bitmap tiltshift_java(Bitmap input, float sigma_far, float sigma_near, int a0, int a1, int a2, int a3){
        Bitmap outBmp = Bitmap.createBitmap(input.getWidth(), input.getHeight(), Bitmap.Config.ARGB_8888);
        //cannot write to input Bitmap, since it may be immutable
        //if you try, you may get a java.lang.IllegalStateException

        int[] pixels = new int[input.getHeight()*input.getWidth()];
        int[] pixelsOut = new int[input.getHeight()*input.getWidth()];
        input.getPixels(pixels,0,input.getWidth(),0,0,input.getWidth(),input.getHeight());

        // Selecting kernel radius as per the sigma values
        int k_radius = (int) ceil(2*max(sigma_far, sigma_near));

        Log.d("SPEEDY_TS", String.valueOf(sigma_far));

        gaussian_filter(pixels,
                pixelsOut,
                0,
                0,
                input.getWidth(),
                a0,
                input.getWidth(),
                input.getHeight(),
                sigma_far,
                k_radius,
                1);

        gaussianGradient_filter(pixels,
                pixelsOut,
                0,
                a0,
                input.getWidth(),
                a1,
                input.getWidth(),
                input.getHeight(),
                sigma_far,
                k_radius,
                0,
                1);

        gaussian_filter(pixels,
                pixelsOut,
                0,
                a1,
                input.getWidth(),
                a2,
                input.getWidth(),
                input.getHeight(),
                0,
                k_radius,
                1);

        gaussianGradient_filter(pixels,
                pixelsOut,
                0,
                a2,
                input.getWidth(),
                a3,
                input.getWidth(),
                input.getHeight(),
                sigma_near,
                k_radius,
                1,
                1);

        gaussian_filter(pixels,
                pixelsOut,
                0,
                a3,
                input.getWidth(),
                input.getHeight(),
                input.getWidth(),
                input.getHeight(),
                sigma_near,
                k_radius,
                1);

        outBmp.setPixels(pixelsOut,0,input.getWidth(),0,0,input.getWidth(),input.getHeight());

        return outBmp;
    }
    public static Bitmap tiltshift_cpp(Bitmap input, float sigma_far, float sigma_near, int a0, int a1, int a2, int a3){
        Bitmap outBmp = Bitmap.createBitmap(input.getWidth(), input.getHeight(), Bitmap.Config.ARGB_8888);
        int[] pixels = new int[input.getHeight()*input.getWidth()];
        int[] pixelsOut = new int[input.getHeight()*input.getWidth()];
        input.getPixels(pixels,0,input.getWidth(),0,0,input.getWidth(),input.getHeight());

        tiltshiftcppnative(pixels,pixelsOut,input.getWidth(),input.getHeight(),sigma_far,sigma_near,a0,a1,a2,a3);

        outBmp.setPixels(pixelsOut,0,input.getWidth(),0,0,input.getWidth(),input.getHeight());
        return outBmp;
    }
    public static Bitmap tiltshift_neon(Bitmap input, float sigma_far, float sigma_near, int a0, int a1, int a2, int a3){
        Bitmap outBmp = Bitmap.createBitmap(input.getWidth(), input.getHeight(), Bitmap.Config.ARGB_8888);
        int[] pixels = new int[input.getHeight()*input.getWidth()];
        int[] pixelsOut = new int[input.getHeight()*input.getWidth()];
        input.getPixels(pixels,0,input.getWidth(),0,0,input.getWidth(),input.getHeight());

        tiltshiftneonnative(pixels,pixelsOut,input.getWidth(),input.getHeight(),sigma_far,sigma_near,a0,a1,a2,a3);

        outBmp.setPixels(pixelsOut,0,input.getWidth(),0,0,input.getWidth(),input.getHeight());
        return outBmp;
    }


    /**
     * A native method that is implemented by the 'native-lib' native library,
     * which is packaged with this application.
     */
    public static native int tiltshiftcppnative(int[] inputPixels, int[] outputPixels, int width, int height, float sigma_far, float sigma_near, int a0, int a1, int a2, int a3);
    public static native int tiltshiftneonnative(int[] inputPixels, int[] outputPixels, int width, int height, float sigma_far, float sigma_near, int a0, int a1, int a2, int a3);

}
