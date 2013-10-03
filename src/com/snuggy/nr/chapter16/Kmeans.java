
package com.snuggy.nr.chapter16;

import static com.snuggy.nr.util.Static.*;

public class Kmeans {

    // Solve for a k-means clustering model from a set of data points and
    // initial guesses of the means. Output is a set of means and an assignment
    // of each data point to one component.
    private int nn, mm, kk, nchg;
    private double[][] data, means;
    private int[] assign, count;

    public Kmeans(final double[][] ddata, final double[][] mmeans) {
        // Constructor. Arguments are the data points (as rows in a matrix),
        // and initial guesses for the means (also as rows in a matrix).
        nn = (nrows(ddata));
        mm = (ncols(ddata));
        kk = (nrows(mmeans));
        data = (ddata);
        means = (mmeans);
        assign = int_arr(nn);
        count = int_arr(kk);
        estep(); // Perform one initial E-step and M-step. User is responsible
        // for calling additional steps until convergence is obtained.
        mstep();
    }
    
    public int[] count() {
        return count;
    }
    
    public int[] assign() {
        return assign;
    }
    
    public double[][] means() {
        return means;
    }

    public int estep() {
        // Perform one E-step.
        int k, m, n, kmin = 0;
        double dmin, d;
        nchg = 0;
        for (k = 0; k < kk; k++)
            count[k] = 0;
        for (n = 0; n < nn; n++) {
            dmin = 9.99e99;
            for (k = 0; k < kk; k++) {
                for (d = 0., m = 0; m < mm; m++)
                    d += SQR(data[n][m] - means[k][m]);
                if (d < dmin) {
                    dmin = d;
                    kmin = k;
                }
            }
            if (kmin != assign[n])
                nchg++;
            assign[n] = kmin;
            count[kmin]++;
        }
        return nchg;
    }

    public void mstep() {
        // Perform one M-step.
        int n, k, m;
        for (k = 0; k < kk; k++)
            for (m = 0; m < mm; m++)
                means[k][m] = 0.;
        for (n = 0; n < nn; n++)
            for (m = 0; m < mm; m++)
                means[assign[n]][m] += data[n][m];
        for (k = 0; k < kk; k++) {
            if (count[k] > 0)
                for (m = 0; m < mm; m++)
                    means[k][m] /= count[k];
        }
    }

}
