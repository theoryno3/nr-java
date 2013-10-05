
package com.snuggy.nr.chapter15;

import static com.snuggy.nr.chapter02.Static.*;
import static com.snuggy.nr.refs.Refs.*;
import static com.snuggy.nr.util.Static.*;

import com.snuggy.nr.refs.*;
import com.snuggy.nr.util.*;

public class Fitlin {

    // Object for general linear least-squares fitting by solving the normal
    // equations, also including the ability to hold specified parameters at
    // fixed, specified values. Call constructor to bind data vectors and
    // fitting functions. Then call any combination of hold, free, and fit as
    // often as desired. fit sets the output quantities a, covar, and chisq.
    private int ndat, ma;
    private final double[] x, y, sig;
    private Func_Doub_To_DoubArr funcs;
    private boolean[] ia;
    private final double[] a; // Output values. a is the vector of fitted
                        // coefficients,
    // covar is its covariance matrix, and chisq is the value of 2 for the fit.
    private final double[][] covar;
    @SuppressWarnings("unused")
    private double chisq;
    
    public final double[] a() {
        return a;
    }
    
    public final double[][] covar() {
        return covar;
    }

    public Fitlin(final double[] xx, final double[] yy, final double[] ssig, final Func_Doub_To_DoubArr funks) {
        // Constructor. Binds references to the data arrays xx, yy, and ssig,
        // and to a user-supplied function funks(x) that returns a VecDoub
        // containing ma basis functions evaluated at x D x. Initializes all
        // parameters as free (not held).
        ndat = (xx.length);
        x = (xx);
        y = (yy);
        sig = (ssig);
        funcs = (funks);
        ma = funcs.eval(x[0]).length;
        a = doub_vec(ma);
        covar = doub_mat(ma, ma);
        ia = bool_vec(ma);
        for (int i = 0; i < ma; i++)
            ia[i] = true;
    }

    public void hold(final int i, final double val) {
        ia[i] = false;
        a[i] = val;
    }

    public void free(final int i) {
        ia[i] = true;
    }

    // Optional functions for holding a parameter, identified by a value i in
    // the range 0; : : : ; ma-1, fixed at the value val, or for freeing a
    // parameter that was previously held fixed. hold and free may be called
    // for any number of parameters before calling fit to calculate best-fit
    // values for the remaining (not held) parameters, and the process may be
    // repeated multiple times. Alternatively, you can set the boolean vector
    // ia directly, before calling fit.

    public void fit() throws NRException {
        // Solve the normal equations for 2 minimization to fit for some or
        // all of the coefficients
        // a[0..ma-1] of a function that depends linearly on a, y D
        // P
        // i ai  funksi .x/. Set answer
        // values for a[0..ma-1], 2 D chisq, and the covariance matrix
        // covar[0..ma-1][0..ma-1]. (Parameters held fixed by calls to hold
        // will return zero covariances.)
        int i, j, k, l, m, mfit = 0;
        double ym, wt, sum, sig2i;
        //final double[] afunc = doub_arr(ma);
        $$double1d afunc = $$(new double[ma]);
        for (j = 0; j < ma; j++)
            if (ia[j])
                mfit++;
        if (mfit == 0)
            throw new NRException("lfit: no parameters to be fitted");
        final double[][] temp = doub_mat(mfit, mfit, 0.), beta = doub_mat(mfit, 1, 0.);
        for (i = 0; i < ndat; i++) { // Loop over data to accumulate
                                     // coefficients of
            $$(afunc, funcs.eval(x[i])); // the normal equations.
            ym = y[i];
            if (mfit < ma) { // Subtract off dependences on known pieces
                for (j = 0; j < ma; j++)
                    // of the fitting function.
                    if (!ia[j])
                        ym -= a[j] * afunc.$()[j];
            }
            sig2i = 1.0 / SQR(sig[i]);
            for (j = 0, l = 0; l < ma; l++) { // Set up matrix and r.h.s. for
                                              // matrix inversion.
                if (ia[l]) {
                    wt = afunc.$()[l] * sig2i;
                    for (k = 0, m = 0; m <= l; m++)
                        if (ia[m])
                            temp[j][k++] += wt * afunc.$()[m];
                    beta[j++][0] += ym * wt;
                }
            }
        }
        for (j = 1; j < mfit; j++)
            for (k = 0; k < j; k++)
                temp[k][j] = temp[j][k];
        gaussj(temp, beta); // Matrix solution.
        for (j = 0, l = 0; l < ma; l++)
            if (ia[l])
                a[l] = beta[j++][0];
        // Spread the solution to appropriate positions in a, and evaluate 2 of
        // the fit.
        chisq = 0.0;
        for (i = 0; i < ndat; i++) {
            $$(afunc, funcs.eval(x[i]));
            sum = 0.0;
            for (j = 0; j < ma; j++)
                sum += a[j] * afunc.$()[j];
            chisq += SQR((y[i] - sum) / sig[i]);
        }
        for (j = 0; j < mfit; j++)
            for (k = 0; k < mfit; k++)
                covar[j][k] = temp[j][k];
        for (i = mfit; i < ma; i++)
            // Rearrange covariance matrix into the correct
            for (j = 0; j < i + 1; j++)
                covar[i][j] = covar[j][i] = 0.0; // order.
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
