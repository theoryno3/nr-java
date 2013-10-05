
package com.snuggy.nr.chapter16;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.chapter02.*;
import com.snuggy.nr.util.*;

public class Gaumixmod extends preGaumixmod {

    // Solve for a Gaussian mixture model from a set of data points and initial
    // guesses of k means.
    private int nn, kk, mm; // Nos. of data points, components, and dimensions.
    private final double[][] data, means, resp; // Local copies of xn’s, k’s, and the
                                          // pnk’s.
    private final double[] frac, lndets; // P.k/’s and log det†k’s.
    private final double[][] sig_arr[]; // †k’s
    private double loglike; // logL.
    
    public final double[][][] sig_arr() {
        return sig_arr;
    }
    
    public final double[][] means() {
        return means;
    }
    
    public final double[] frac() {
        return frac;
    }

    public Gaumixmod(final double[][] ddata, final double[][] mmeans) throws NRException {
        // Constructor. Arguments are the data points (as rows in a matrix)
        // and initial guesses for the means (also as rows in a matrix).
        super(ncols(ddata));
        nn = (nrows(ddata));
        kk = (nrows(mmeans));
        mm = (mmstat);
        data = (ddata);
        means = (mmeans);
        resp = doub_mat(nn, kk);
        frac = doub_vec(kk);
        lndets = doub_vec(kk);
        sig_arr = new double[kk][][];
        int i, j, k;
        for (k = 0; k < kk; k++) {
            sig_arr[k] = Mat_mm();
            frac[k] = 1. / kk; // Uniform prior on P.k/.
            for (i = 0; i < mm; i++) {
                for (j = 0; j < mm; j++)
                    sig_arr[k][i][j] = 0.;
                sig_arr[k][i][i] = 1.0e-10; // See text at end of this section.
            }
        }
        estep(); // Perform one initial E-step and M-step. User
        // is responsible for calling additional steps until convergence is
        // obtained.
        mstep();
    }

    public double estep() throws NRException {
        // Perform one E-step of the EM algorithm.
        int k, m, n;
        double tmp, sum, max, oldloglike;
        final double[] u = doub_vec(mm), v = doub_vec(mm);
        oldloglike = loglike;
        for (k = 0; k < kk; k++) { // Outer loop for computing the pnk’s.
            Cholesky choltmp = new Cholesky(sig_arr[k]); // Decompose †k in the
                                                     // outer loop.
            lndets[k] = choltmp.logdet();
            for (n = 0; n < nn; n++) { // Inner loop for pnk’s.
                for (m = 0; m < mm; m++)
                    u[m] = data[n][m] - means[k][m];
                choltmp.elsolve(u, v); // Solve L  v D u.
                for (sum = 0., m = 0; m < mm; m++)
                    sum += SQR(v[m]);
                resp[n][k] = -0.5 * (sum + lndets[k]) + log(frac[k]);
            }
        }
        // At this point we have unnormalized logs of the pnk’s. We need to
        // normalize using log-sum-exp and compute the log-likelihood.
        loglike = 0;
        for (n = 0; n < nn; n++) { // Separate normalization for each n.
            max = -99.9e99; // Log-sum-exp trick begins here.
            for (k = 0; k < kk; k++)
                if (resp[n][k] > max)
                    max = resp[n][k];
            for (sum = 0., k = 0; k < kk; k++)
                sum += exp(resp[n][k] - max);
            tmp = max + log(sum);
            for (k = 0; k < kk; k++)
                resp[n][k] = exp(resp[n][k] - tmp);
            loglike += tmp;
        }
        return loglike - oldloglike; // When abs of this is small, then we have
    } // converged.

    public void mstep() {
        // Perform one M-step of the EM algorithm.
        int j, n, k, m;
        double wgt, sum;
        for (k = 0; k < kk; k++) {
            wgt = 0.;
            for (n = 0; n < nn; n++)
                wgt += resp[n][k];
            frac[k] = wgt / nn; // Equation (16.1.7).
            for (m = 0; m < mm; m++) {
                for (sum = 0., n = 0; n < nn; n++)
                    sum += resp[n][k] * data[n][m];
                means[k][m] = sum / wgt; // Equation (16.1.6).
                for (j = 0; j < mm; j++) {
                    for (sum = 0., n = 0; n < nn; n++) {
                        sum += resp[n][k] * (data[n][m] - means[k][m]) * (data[n][j] - means[k][j]);
                    }
                    sig_arr[k][m][j] = sum / wgt; // Equation (16.1.6).
                }
            }
        }
    }

}
