
package com.snuggy.nr.chapter15;

import static com.snuggy.nr.refs.Refs.*;
import static com.snuggy.nr.util.Static.*;

import com.snuggy.nr.chapter02.*;
import com.snuggy.nr.refs.*;
import com.snuggy.nr.util.*;

@Deprecated @Broken
public class Fitsvd {

    // Object for general linear least-squares fitting using singular value
    // decomposition. Call constructor to bind data vectors and fitting
    // functions. Then call fit, which sets the output quantities a, covar,
    // and chisq.
    private int ndat, ma;
    private double tol;
    private double[] x_arr[]; // (Why is x a pointer? Explained in 15.4.4.)
    private int x_off;
    private double[] y, sig; // (Why is x a pointer? Explained in 15.4.4.)
    private Func_Doub_To_DoubArr funcs;
    private double[] a; // Output values. a is the vector of fitted
                        // coefficients,
    // covar is its covariance matrix, and chisq is the value of 2 for the fit.
    private double[][] covar;
    @SuppressWarnings("unused")
    private double chisq;
    
    public double[][] covar() {
        return covar;
    }
    
    public double[] a() {
        return a;
    }

    public Fitsvd(final double[] xx_arr[], final int xx_off, final double[] yy, final double[] ssig, Func_Doub_To_DoubArr funks) {
        this(xx_arr, xx_off, yy, ssig, funks, 1.e-12);
    }

    public Fitsvd(final double[] xx_arr[], final int xx_off, final double[] yy, final double[] ssig, Func_Doub_To_DoubArr funks,
            final double TOL) {
        // Constructor. Binds references to the data arrays xx, yy, and ssig,
        // and to a user-supplied function funks(x) that returns a VecDoub
        // containing ma basis functions evaluated at x D x. If TOL is positive,
        // it is the threshold (relative to the largest singular value) for
        // discarding small singular values. If it is  0, the default value in
        // SVD is used.
        ndat = (yy.length);
        // x = (&xx);
        x_arr = xx_arr;
        x_off = xx_off;
        xmd_arr = null;
        xmd_off = 0;
        y = (yy);
        sig = (ssig);
        funcs = (funks);
        tol = (TOL);
    }

    public void fit() throws NRException {
        // Solve by singular value decomposition the 2 minimization that fits
        // for the coefficients a[0..ma-1] of a function that depends linearly
        // on a, y D
        // P
        // i ai  funksi .x/. Set answer values for a[0..ma-1], chisq D 2, and
        // the covariance matrix covar[0..ma-1][0..ma-1].
        int i, j, k;
        double tmp, thresh, sum;
        if (x_arr != null)
            ma = funcs.eval(x_arr[x_off][0]).length;
        else
            ma = funcsmd.eval(row(xmd_arr[xmd_off], 0)).length; // (Discussed in
                                                                // 15.4.4.)
        a = doub_arr(ma);
        covar = doub_mat(ma, ma);
        double[][] aa = doub_mat(ndat, ma);
        double[] b = doub_arr(ndat);
        $double1d afunc = $(new double[ma]); 
        for (i = 0; i < ndat; i++) { // Accumulate coefficients of the
            if (x_arr != null)
                $$(afunc, funcs.eval((x_arr[x_off])[i])); // design matrix.
            else
                $$(afunc, funcsmd.eval(row(xmd_arr[xmd_off], i))); // (Discussed in
                                                                // 15.4.4.)
            tmp = 1.0 / sig[i];
            for (j = 0; j < ma; j++)
                aa[i][j] = afunc.$()[j] * tmp;
            b[i] = y[i] * tmp;
        }
        SVD svd = new SVD(aa); // Singular value decomposition.
        thresh = (tol > 0. ? tol * svd.w()[0] : -1.);
        svd.solve(b, a, thresh); // Solve for the coefficients.
        chisq = 0.0; // Evaluate chi-square.
        for (i = 0; i < ndat; i++) {
            sum = 0.;
            for (j = 0; j < ma; j++)
                sum += aa[i][j] * a[j];
            chisq += SQR(sum - b[i]);
        }
        for (i = 0; i < ma; i++) { // Sum contributions to covariance
            for (j = 0; j < i + 1; j++) { // matrix (15.4.20).
                sum = 0.0;
                for (k = 0; k < ma; k++)
                    if (svd.w()[k] > svd.tsh())
                        sum += svd.v()[i][k] * svd.v()[j][k] / SQR(svd.w()[k]);
                covar[j][i] = covar[i][j] = sum;
            }
        }
    }

    // From here on, code for multidimensional fits, to be discussed in 15.4.4.

    private double[][] xmd_arr[];
    private int xmd_off;

    // VecDoub (*funcsmd)(VecDoub_I &);
    private Func_DoubArr_To_DoubArr funcsmd;

    public Fitsvd(final double[][] xx_arr[], final int xx_off, final double[] yy, final double[] ssig, Func_DoubArr_To_DoubArr funks) {
        this(xx_arr, xx_off, yy, ssig, funks, 1.e-12);
    }

    public Fitsvd(final double[][] xx_arr[], final int xx_off, final double[] yy, final double[] ssig, Func_DoubArr_To_DoubArr funks,
            final double TOL) {
        // Constructor for multidimensional fits. Exactly the same as the
        // previous constructor, except that xx is now a matrix whose rows are
        // the multidimensional data points and funks is now a function of a
        // multidimensional data point (as a VecDoub).
        ndat = (yy.length);
        x_arr = (null);
        x_off = 0;
        xmd_arr = (xx_arr);
        xmd_off = xx_off;
        y = (yy);
        sig = (ssig);
        funcsmd = (funks);
        tol = (TOL);
    }

    private double[] row(final double[][] a, final int i) {
        // Utility. Returns the row of a MatDoub as a VecDoub.
        int j, n = ncols(a);
        double[] ans = doub_arr(n);
        for (j = 0; j < n; j++)
            ans[j] = a[i][j];
        return ans;
    }

}
