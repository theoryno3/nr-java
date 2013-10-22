
package com.snuggy.nr.chapter15;

import static com.snuggy.nr.chapter02.Static.*;
import static com.snuggy.nr.refs.Refs.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.refs.*;
import com.snuggy.nr.util.*;

public class Fitmrq {

    // Object for nonlinear least-squares fitting by the Levenberg-Marquardt
    // method, also including the ability to hold specified parameters at
    // fixed, specified values. Call constructor to bind data vectors and
    // fitting functions and to input an initial parameter guess. Then call
    // any combination of hold, free, and fit as often as desired. fit sets
    // the output quantities a, covar, alpha, and chisq.
    private static final int NDONE = 4, ITMAX = 1000; // Convergence parameters.

    private int ndat, ma, mfit;
    private final double[] x, y, sig;
    private double tol;
    private Func_Doub_DoubVec_DoubRef_DoubVec_To_Void funcs;
    private boolean[] ia;
    private final double[] a; // Output values. a is the vector of fitted
                        // coefficients,
    // covar is its covariance matrix, alpha is the curvature matrix, and
    // chisq is the value of 2 for the fit.
    private final double[][] covar;
    private final double[][] alpha;
    private double chisq;
    
    public final double[] a() {
        return a;
    }
    
    public final double[][] covar() {
        return covar;
    }

    public Fitmrq(final double[] xx, final double[] yy, final double[] ssig, final double[] aa, Func_Doub_DoubVec_DoubRef_DoubVec_To_Void funks) {
        this(xx, yy, ssig, aa, funks, 1.e-3);
    }

    public Fitmrq(final double[] xx, final double[] yy, final double[] ssig, final double[] aa, Func_Doub_DoubVec_DoubRef_DoubVec_To_Void funks,
            final double TOL) {
        // Constructor. Binds references to the data arrays xx, yy, and ssig,
        // and to a user-supplied function funks that calculates the nonlinear
        // fitting function and its derivatives. Also inputs the initial
        // parameters guess aa (which is copied, not modified) and an optional
        // convergence tolerance TOL. Initializes all parameters as free (not
        // held).
        ndat = (xx.length);
        ma = (aa.length);
        x = (xx);
        y = (yy);
        sig = (ssig);
        tol = (TOL);
        funcs = (funks);
        ia = bool_vec(ma);
        alpha = doub_mat(ma, ma);
        a = (aa);
        covar = doub_mat(ma, ma);
        for (int i = 0; i < ma; i++)
            ia[i] = true;
    }

    public void hold(final int i, final double val) {
        ia[i] = false;
        a[i] = val;
    }

    public void free(final int i) {
        // Optional functions for holding a parameter, identified by a value i
        // in the range 0; : : : ;ma-1, fixed at the value val, or for freeing
        // a parameter that was previously held fixed. hold and free may be
        // called for any number of parameters before calling fit to calculate
        // best-fit values for the remaining (not held) parameters, and the
        // process may be repeated multiple times.
        ia[i] = true;
    }

    public void fit() throws NRException {
        // Iterate to reduce the 2 of a fit between a set of data points
        // x[0..ndat-1], y[0..ndat-1] with individual standard deviations
        // sig[0..ndat-1], and a nonlinear function that depends on ma
        // coefficients a[0..ma-1]. When 2 is no longer decreasing, set
        // best-fit values for the parameters a[0..ma-1], and chisq D 2,
        // covar[0..ma-1][0..ma-1], and alpha[0..ma-1][0..ma-1].
        // (Parameters held fixed will return zero covariances.)
        int j, k, l, iter, done = 0;
        double alamda = .001, ochisq;
        final double[] atry = doub_vec(ma), beta = doub_vec(ma), da = doub_vec(ma);
        mfit = 0;
        for (j = 0; j < ma; j++)
            if (ia[j])
                mfit++;
        final double[][] oneda = doub_mat(mfit, 1), temp = doub_mat(mfit, mfit);
        mrqcof(a, alpha, beta); // Initialization.
        for (j = 0; j < ma; j++)
            atry[j] = a[j];
        ochisq = chisq;
        for (iter = 0; iter < ITMAX; iter++) {
            if (done == NDONE)
                alamda = 0.; // Last pass. Use zero alamda.
            for (j = 0; j < mfit; j++) { // Alter linearized fitting matrix, by
                                         // augmenting di
                for (k = 0; k < mfit; k++)
                    covar[j][k] = alpha[j][k]; // agonal elements.
                covar[j][j] = alpha[j][j] * (1.0 + alamda);
                for (k = 0; k < mfit; k++)
                    temp[j][k] = covar[j][k];
                oneda[j][0] = beta[j];
            }
            gaussj(temp, oneda); // Matrix solution.
            for (j = 0; j < mfit; j++) {
                for (k = 0; k < mfit; k++)
                    covar[j][k] = temp[j][k];
                da[j] = oneda[j][0];
            }
            if (done == NDONE) { // Converged. Clean up and return.
                covsrt(covar);
                covsrt(alpha);
                return;
            }
            for (j = 0, l = 0; l < ma; l++)
                // Did the trial succeed?
                if (ia[l])
                    atry[l] = a[l] + da[j++];
            mrqcof(atry, covar, da);
            if (abs(chisq - ochisq) < MAX(tol, tol * chisq))
                done++;
            if (chisq < ochisq) { // Success, accept the new solution.
                alamda *= 0.1;
                ochisq = chisq;
                for (j = 0; j < mfit; j++) {
                    for (k = 0; k < mfit; k++)
                        alpha[j][k] = covar[j][k];
                    beta[j] = da[j];
                }
                for (l = 0; l < ma; l++)
                    a[l] = atry[l];
            } else { // Failure, increase alamda.
                alamda *= 10.0;
                chisq = ochisq;
            }
        }
        throw new NRException("Fitmrq too many iterations");
    }

    public void mrqcof(final double[] a, final double[][] alpha, final double[] beta) {
        // Used by fit to evaluate the linearized fitting matrix alpha, and
        // vector beta as in (15.5.8), and to calculate 2.
        int i, j, k, l, m;
        double wt, sig2i, dy;
        final $double ymod = $(0.0);
        final double[] dyda = doub_vec(ma);
        for (j = 0; j < mfit; j++) { // Initialize (symmetric) alpha, beta.
            for (k = 0; k <= j; k++)
                alpha[j][k] = 0.0;
            beta[j] = 0.;
        }
        chisq = 0.;
        for (i = 0; i < ndat; i++) { // Summation loop over all data.
            funcs.eval(x[i], a, ymod, dyda);
            sig2i = 1.0 / (sig[i] * sig[i]);
            dy = y[i] - ymod.$();
            for (j = 0, l = 0; l < ma; l++) {
                if (ia[l]) {
                    wt = dyda[l] * sig2i;
                    for (k = 0, m = 0; m < l + 1; m++)
                        if (ia[m])
                            alpha[j][k++] += wt * dyda[m];
                    beta[j++] += dy * wt;
                }
            }
            chisq += dy * dy * sig2i; // And find 2.
        }
        for (j = 1; j < mfit; j++)
            // Fill in the symmetric side.
            for (k = 0; k < j; k++)
                alpha[k][j] = alpha[j][k];
    }

    public void covsrt(final double[][] covar) {
        // Expand in storage the covariance matrix covar, so as to take into
        // account parameters that are being held fixed. (For the latter,
        // return zero covariances.)
        int i, j, k;
        for (i = mfit; i < ma; i++)
            for (j = 0; j < i + 1; j++)
                covar[i][j] = covar[j][i] = 0.0;
        k = mfit - 1;
        for (j = ma - 1; j >= 0; j--) {
            if (ia[j]) {
                for (i = 0; i < ma; i++)
                    SWAP(covar, i, k, i, j);
                for (i = 0; i < ma; i++)
                    SWAP(covar, k, i, j, i);
                k--;
            }
        }
    }

}
