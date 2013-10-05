
package com.snuggy.nr.chapter03;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.chapter02.*;
import com.snuggy.nr.util.*;

public class RBF_interp {

    // Object for radial basis function interpolation using n points in dim
    // dimensions. Call constructor once, then interp as many times as desired.

    private int dim, n;
    private final double[][] pts;
    private final double[] vals;
    private final double[] w;
    private RBF_fn fn;
    private boolean norm;

    public RBF_interp(final double[][] ptss, final double[] valss, final RBF_fn func) 
            throws NRException {
        this(ptss, valss, func, false);
    }

    public RBF_interp(final double[][] ptss, final double[] valss, final RBF_fn func, final boolean nrbf) 
            throws NRException {
        // Constructor. The n  dim matrix ptss inputs the data points, the
        // vector valss the function values. func contains the chosen radial
        // basis function, derived from the class RBF_fn. The default value
        // of nrbf gives RBF interpolation; set it to 1 for NRBF.
        dim = (ncols(ptss));
        n = (nrows(ptss));
        pts = (ptss);
        vals = (valss);
        w = doub_arr(n);
        fn = (func);
        norm = (nrbf);
        int i, j;
        double sum;
        final double[][] rbf = doub_mat(n, n);
        final double[] rhs = doub_arr(n);
        for (i = 0; i < n; i++) { // Fill the matrix .jri rj j/ and the r.h.s.
                                  // vector.
            sum = 0.;
            for (j = 0; j < n; j++) {
                sum += (rbf[i][j] = fn.rbf(rad(pts[i], 0, pts[j], 0)));
            }
            if (norm)
                rhs[i] = sum * vals[i];
            else
                rhs[i] = vals[i];
        }
        LUdcmp lu = new LUdcmp(rbf); // Solve the set of linear equations.
        lu.solve(rhs, w);
    }

    public double interp(final double[] pt) throws NRException {
        // Return the interpolated function value at a dim-dimensional point pt.
        double fval, sum = 0., sumw = 0.;
        if (pt.length != dim)
            throw new NRException("RBF_interp bad pt size");
        for (int i = 0; i < n; i++) { // Sum over all tabulated points.
            fval = fn.rbf(rad(pt, 0, pts[i], 0));
            sumw += w[i] * fval;
            sum += fval;
        }
        return norm ? sumw / sum : sumw;
    }

    public double rad(final double[] p1_arr, final int p1_off, 
                        final double[] p2_arr, final int p2_off) {
        // Euclidean distance.
        double sum = 0.;
        for (int i = 0; i < dim; i++)
            sum += SQR(p1_arr[p1_off + i] - p2_arr[p2_off + i]);
        return sqrt(sum);
    }

}
