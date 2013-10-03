
package com.snuggy.nr.chapter16;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

public class Svmpolykernel extends Svmgenkernel {

    // Kernel structure for the polynomial kernel.
    private int n;
    private double a, b, d;

    public Svmpolykernel(final double[][] ddata, final double[] yy, final double aa, final double bb, final double dd) {
        // Constructor is called with the m  n data matrix, the vector of yi
        // ’s,
        // length m, and the constants a, b, and d.
        super(yy, ddata);
        n = (ncols(data));
        a = (aa);
        b = (bb);
        d = (dd);
        fill();
    }

    public double kernel(final double xi_arr[], final int xi_off, final double xj_arr[], final int xj_off) {
        double dott = 0.;
        for (int k = 0; k < n; k++)
            dott += xi_arr[xi_off + k] * xj_arr[xj_off + k];
        return pow(a * dott + b, d);
    }

}
