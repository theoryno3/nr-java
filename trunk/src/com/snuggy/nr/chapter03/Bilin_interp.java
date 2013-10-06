
package com.snuggy.nr.chapter03;

import static com.snuggy.nr.util.Static.*;

import com.snuggy.nr.util.*;

public class Bilin_interp {

    // Object for bilinear interpolation on a matrix. Construct with a
    // vector of x1 values, a vector of x2 values, and a matrix of tabulated
    // function values yij. Then call interp for interpolated values.

    protected int m, n;
    protected final double[][] y;
    protected Linear_interp x1terp, x2terp;

    public Bilin_interp(final double[] x1v, final double[] x2v, final double[][] ym) {
        m = (x1v.length);
        n = (x2v.length);
        y = doub_mat(ym);
        x1terp = new Linear_interp(x1v, x1v);
        x2terp = new Linear_interp(x2v, x2v);
    } // Construct dummy 1-dim interpolations
      // for their locate and hunt

    public double interp(final double x1p, final double x2p) throws NRException { // functions.
        int i, j;
        double yy, t, u;
        i = (x1terp.cor != 0) ? x1terp.hunt(x1p) : x1terp.locate(x1p);
        j = (x2terp.cor != 0) ? x2terp.hunt(x2p) : x2terp.locate(x2p);
        // Find the grid square.
        t = (x1p - x1terp.xx.$_(i))
                / (x1terp.xx.$_(i + 1) - x1terp.xx.$_(i)); // Interpolate.
        u = (x2p - x2terp.xx.$_(j))
                / (x2terp.xx.$_(j + 1) - x2terp.xx.$_(j));
        yy = (1. - t) * (1. - u) * y[i][j] + t * (1. - u) * y[i + 1][j] + (1. - t) * u * y[i][j + 1] + t * u
                * y[i + 1][j + 1];
        return yy;
    }

}
