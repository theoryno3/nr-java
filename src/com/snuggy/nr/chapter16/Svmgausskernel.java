
package com.snuggy.nr.chapter16;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

public class Svmgausskernel extends Svmgenkernel {

    // Kernel structure for the Gaussian radial basis function kernel.
    private int n;
    private double sigma;

    public Svmgausskernel(final double[][] ddata, final double[] yy, final double ssigma) {
        // Constructor is called with the m  n data matrix, the vector of
        // yi ’s, length m, and the constant .
        super(yy, ddata);
        n = (ncols(data));
        sigma = (ssigma);
        fill();
    }

    public double kernel(final double xi_arr[], final int xi_off, final double xj_arr[], final int xj_off) {
        double dott = 0.;
        for (int k = 0; k < n; k++)
            dott += SQR(xi_arr[xi_off + k] - xj_arr[xj_off + k]);
        return exp(-0.5 * dott / (sigma * sigma));
    }

}
