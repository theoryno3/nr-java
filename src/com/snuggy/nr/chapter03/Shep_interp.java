
package com.snuggy.nr.chapter03;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Shep_interp {

    // Object for Shepard interpolation using n points in dim dimensions.
    // Call constructor once, then interp as many times as desired.

    private int dim, n;
    private final double[][] pts;
    private final double[] vals;
    private double pneg;

    public Shep_interp(final double[][] ptss, final double[] valss) {
        this(ptss, valss, 2.0);
    }

    public Shep_interp(final double[][] ptss, final double[] valss, final double p) {
        dim = (ncols(ptss));
        n = (nrows(ptss));
        pts = (ptss);
        vals = (valss);
        pneg = (-p);
    }

    // Constructor. The n  dim matrix ptss inputs the data points, the
    // vector valss the function values. Set p to the desired exponent.
    // The default value is typical.

    public double interp(final double[] pt) throws NRException {
        // Return the interpolated function value at a dim-dimensional point pt.
        double r, w, sum = 0., sumw = 0.;
        if (pt.length != dim)
            throw new NRException("RBF_interp bad pt size");
        for (int i = 0; i < n; i++) {
            if ((r = rad(pt, 0, pts[i], 0)) == 0.)
                return vals[i];
            sum += (w = pow(r, pneg));
            sumw += w * vals[i];
        }
        return sumw / sum;
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
