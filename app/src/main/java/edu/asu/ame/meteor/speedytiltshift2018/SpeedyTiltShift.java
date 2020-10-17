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

    public static Bitmap tiltshift_java(Bitmap input, float sigma_far, float sigma_near, int a0, int a1, int a2, int a3){
        Bitmap outBmp = Bitmap.createBitmap(input.getWidth(), input.getHeight(), Bitmap.Config.ARGB_8888);
        //cannot write to input Bitmap, since it may be immutable
        //if you try, you may get a java.lang.IllegalStateException

        int[] pixels = new int[input.getHeight()*input.getWidth()];
        int[] pixelsOut = new int[input.getHeight()*input.getWidth()];
        input.getPixels(pixels,0,input.getWidth(),0,0,input.getWidth(),input.getHeight());

        // Selecting kernel radius as per the sigma values
        int k_radius = (int) ceil(2*max(sigma_far, sigma_near));
        int fast = 1;

//        sigma_far = 5.0f;
//        sigma_near = 5.0f;

//        a0 = input.getHeight()/5;
//        a1 = (input.getHeight()/5)*2;
//        a1 = (input.getHeight()/5)*3;
//        a1 = (input.getHeight()/5)*4;

        // Log.d("SPEEDY_TS", String.valueOf(a0) + ", " + String.valueOf(a1) + ", " + String.valueOf(a2) + ", " + String.valueOf(a3));

        SpeedyTiltShiftJava strip0 = new SpeedyTiltShiftJava(pixels,
                pixelsOut,
                0,
                0,
                input.getWidth(),
                a0,
                input.getWidth(),
                input.getHeight(),
                sigma_far,
                k_radius,
                2,
                fast);

        SpeedyTiltShiftJava strip1 = new SpeedyTiltShiftJava(pixels,
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
                fast);

        SpeedyTiltShiftJava strip2 = new SpeedyTiltShiftJava(pixels,
                pixelsOut,
                0,
                a1,
                input.getWidth(),
                a2,
                input.getWidth(),
                input.getHeight(),
                0,
                k_radius,
                2,
                fast);

        SpeedyTiltShiftJava strip3 = new SpeedyTiltShiftJava(pixels,
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
                fast);

        SpeedyTiltShiftJava strip4 = new SpeedyTiltShiftJava(pixels,
                pixelsOut,
                0,
                a3,
                input.getWidth(),
                input.getHeight(),
                input.getWidth(),
                input.getHeight(),
                sigma_near,
                k_radius,
                2,
                fast);

//        strip0.printInfo(0);
//        strip1.printInfo(1);
//        strip2.printInfo(2);
//        strip3.printInfo(3);
//        strip4.printInfo(4);

        strip0.start();
        strip1.start();
        strip2.start();
        strip3.start();
        strip4.start();

        try {
            strip0.join();
        } catch (InterruptedException e) {
            Log.d("SPEEDY_TS", "Error !");
        }
        try {
            strip1.join();
        } catch (InterruptedException e) {
            Log.d("SPEEDY_TS", "Error !");
        }
        try {
            strip2.join();
        } catch (InterruptedException e) {
            Log.d("SPEEDY_TS", "Error !");
        }
        try {
            strip3.join();
        } catch (InterruptedException e) {
            Log.d("SPEEDY_TS", "Error !");
        }
        try {
            strip4.join();
        } catch (InterruptedException e) {
            Log.d("SPEEDY_TS", "Error !");
        }

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
