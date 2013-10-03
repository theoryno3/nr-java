
package com.snuggy.nr.chapter03;

import static com.snuggy.nr.util.Static.*;

import com.snuggy.nr.util.*;

public class Poly2D_interp {

    // Object for two-dimensional polynomial interpolation on a matrix.
    // Construct with a vector of x1 values, a vector of x2 values, a
    // matrix of tabulated function values yij , and integers to specify the
    // number of points to use locally in each direction. Then call interp
    // for interpolated values.

    protected int m, n, mm, nn;
    protected final double[][] y;
    protected double[] yv;
    protected Poly_interp x1terp, x2terp;

    public Poly2D_interp(final double[] x1v, final double[] x2v, final double[][] ym, final int mp, final int np) {
        m = (x1v.length);
        n = (x2v.length);
        mm = (mp);
        nn = (np);
        y = doub_mat(ym);
        yv = doub_arr(m);
        x1terp = new Poly_interp(x1v, yv, mm);
        x2terp = new Poly_interp(x2v, x2v, nn);
    }

    // Dummy 1-dim interpolations for their
    // locate and hunt functions.

    public double interp(final double x1p, final double x2p) throws NRException {
        int i, j, k;
        i = (x1terp.cor != 0) ? x1terp.hunt(x1p) : x1terp.locate(x1p);
        j = (x2terp.cor != 0) ? x2terp.hunt(x2p) : x2terp.locate(x2p);
        // Find grid block.
        for (k = i; k < i + mm; k++) { // mm interpolations in the x2 direction.
            // x2terp.yy = &y[k][0];
            x2terp.yy_arr = y[k];
            x2terp.yy_off = 0;
            yv[k] = x2terp.rawinterp(j, x2p);
        }
        return x1terp.rawinterp(i, x1p);
    } // A final interpolation in the x1 direc
      // tion.

}
