
package com.snuggy.nr.chapter16;

import static com.snuggy.nr.util.Static.*;

public abstract class Svmgenkernel {

    // Virtual class that defines what a kernel structure needs to provide.
    protected int m; // No. of data points; counter for kernel calls.
    @SuppressWarnings("unused")
    private int kcalls;
    private final double[][] ker; // Locally stored kernel matrix.
    private final double[] y; // Must provide reference to the yi ’s.
    protected final double[][] data; // Must provide reference to the xi ’s.

    public Svmgenkernel(final double[] yy, final double[][] ddata) {
        // Every kernel structure must provide a kernel function that returns
        // the kernel for arbitrary feature vectors.
        m = (yy.length);
        kcalls = (0);
        ker = doub_mat(m, m);
        y = (yy);
        data = (ddata);
    }
    
    public final double[][] ker() {
        return ker;
    }
    
    public abstract double kernel(final double xi_arr[], final int xi_off, final double xj_arr[], final int xj_off);

    public double kernel(final int i, final double xj_arr[], final int xj_off) {
        return kernel(data[i], 0, xj_arr, xj_off);
    }

    // Every kernel structure’s constructor must call fill to fill the ker
    // matrix.

    void fill() {
        int i, j;
        for (i = 0; i < m; i++)
            for (j = 0; j <= i; j++) {
                ker[i][j] = ker[j][i] = kernel(data[i], 0, data[j], 0);
            }
    }

    public final double[] y() {
        return y;
    }

    /*
    public void set_t(final double[] y) {
        this.y = y;
    }
    */

}
