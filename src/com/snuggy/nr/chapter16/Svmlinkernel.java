
package com.snuggy.nr.chapter16;

import static com.snuggy.nr.util.Static.*;

public class Svmlinkernel extends Svmgenkernel {

    // Kernel structure for the linear kernel, the dot product of two feature
    // vectors (with overall means of each component subtracted).
    private int n;
    private final double[] mu;

    public Svmlinkernel(final double[][] ddata, final double[] yy) {
        // Constructor is called with the m n data matrix, and the vector
        // of yi ’s, length m.
        super(yy, ddata);
        n = (ncols(data));
        mu = doub_vec(n);
        int i, j;
        for (j = 0; j < n; j++)
            mu[j] = 0.;
        for (i = 0; i < m; i++)
            for (j = 0; j < n; j++)
                mu[j] += data[i][j];
        for (j = 0; j < n; j++)
            mu[j] /= m;
        fill();
    }

    public double kernel(final double xi_arr[], final int xi_off, final double xj_arr[], final int xj_off) {
        double dott = 0.;
        for (int k = 0; k < n; k++)
            dott += (xi_arr[xi_off + k] - mu[k]) * (xj_arr[xj_off + k] - mu[k]);
        return dott;
    }

}
