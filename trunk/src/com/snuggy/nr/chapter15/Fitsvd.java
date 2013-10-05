
package com.snuggy.nr.chapter15;

import static com.snuggy.nr.refs.Refs.*;
import static com.snuggy.nr.util.Static.*;

import com.snuggy.nr.chapter02.*;
import com.snuggy.nr.refs.*;
import com.snuggy.nr.util.*;

public class Fitsvd {

    // Object for general linear least-squares fitting using singular value
    // decomposition. Call constructor to bind data vectors and fitting
    // functions. Then call fit, which sets the output quantities a, covar,
    // and chisq.
    private int ndat, ma;
    private double tol;
    private double[] x; // (Why is x a pointer? Explained in 15.4.4.)
    private final double[] y, sig; 
    private Func_Doub_To_DoubArr funcs;
    private $$double1d a; // Output values. a is the vector of fitted
                        // coefficients,
    // covar is its covariance matrix, and chisq is the value of 2 for the fit.
    private $$double2d covar;
    @SuppressWarnings("unused")
    private double chisq;
    
    public Fitsvd(final double[] xx, final double[] yy, final double[] ssig, Func_Doub_To_DoubArr funks) throws NRException {
        this(xx, yy, ssig, funks, 1.e-12);
    }
    
    public double[] a() {
        return a.$();
    }
    
    public double[][] covar() {
        return covar.$();
    }

    public Fitsvd(final double[] xx, final double[] yy, final double[] ssig, Func_Doub_To_DoubArr funks,
            final double TOL) throws NRException {
        // Constructor. Binds references to the data arrays xx, yy, and ssig,
        // and to a user-supplied function funks(x) that returns a VecDoub
        // containing ma basis functions evaluated at x D x. If TOL is positive,
        // it is the threshold (relative to the largest singular value) for
        // discarding small singular values. If it is  0, the default value in
        // SVD is used.
        ndat = (yy.length);
        // x = (&xx);
        x = xx;
        xmd = null;
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
        if (x != null)
            ma = funcs.eval(x[0]).length;
        else
            ma = funcsmd.eval(row(xmd.$(), 0)).length; // (Discussed in
                                                                // 15.4.4.)
        a = $$(doub_arr(ma));
        covar = $$(doub_mat(ma, ma));
        final double[][] aa = doub_mat(ndat, ma);
        final double[] b = doub_arr(ndat);
        $$double1d afunc = $$(new double[ma]); 
        for (i = 0; i < ndat; i++) { // Accumulate coefficients of the
            if (x != null)
                $$(afunc, funcs.eval(x[i])); // design matrix.
            else
                $$(afunc, funcsmd.eval(row(xmd.$(), i))); // (Discussed in
                                                                // 15.4.4.)
            tmp = 1.0 / sig[i];
            for (j = 0; j < ma; j++)
                aa[i][j] = afunc.$()[j] * tmp;
            b[i] = y[i] * tmp;
        }
        SVD svd = new SVD(aa); // Singular value decomposition.
        thresh = (tol > 0. ? tol * svd.w()[0] : -1.);
        svd.solve(b, a.$(), thresh); // Solve for the coefficients.
        chisq = 0.0; // Evaluate chi-square.
        for (i = 0; i < ndat; i++) {
            sum = 0.;
            for (j = 0; j < ma; j++)
                sum += aa[i][j] * a.$()[j];
            chisq += SQR(sum - b[i]);
        }
        for (i = 0; i < ma; i++) { // Sum contributions to covariance
            for (j = 0; j < i + 1; j++) { // matrix (15.4.20).
                sum = 0.0;
                for (k = 0; k < ma; k++)
                    if (svd.w()[k] > svd.tsh())
                        sum += svd.v()[i][k] * svd.v()[j][k] / SQR(svd.w()[k]);
                $(covar, j, i, $(covar, i, j, sum));
            }
        }
    }

    // From here on, code for multidimensional fits, to be discussed in 15.4.4.

    private $$double2d xmd;

    // VecDoub (*funcsmd)(VecDoub_I &);
    private Func_DoubArr_To_DoubArr funcsmd;

    public Fitsvd(final $$double2d xx, final double[] yy, final double[] ssig, Func_DoubArr_To_DoubArr funks) throws NRException {
        this(xx, yy, ssig, funks, 1.e-12);
    }

    public Fitsvd(final $$double2d xx, final double[] yy, final double[] ssig, Func_DoubArr_To_DoubArr funks,
            final double TOL) throws NRException {
        // Constructor for multidimensional fits. Exactly the same as the
        // previous constructor, except that xx is now a matrix whose rows are
        // the multidimensional data points and funks is now a function of a
        // multidimensional data point (as a VecDoub).
        ndat = (yy.length);
        x = (null);
        xmd.$(xx.$());
        y = (yy);
        sig = (ssig);
        funcsmd = (funks);
        tol = (TOL);
    }

    private final double[] row(final double[][] a, final int i) {
        // Utility. Returns the row of a MatDoub as a VecDoub.
        int j, n = ncols(a);
        final double[] ans = doub_arr(n);
        for (j = 0; j < n; j++)
            ans[j] = a[i][j];
        return ans;
    }

}
