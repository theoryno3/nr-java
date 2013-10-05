
package com.snuggy.nr.chapter03;

import static com.snuggy.nr.util.Static.*;

import java.util.*;

import com.snuggy.nr.util.*;

public class Spline2D_interp {

    // Object for two-dimensional cubic spline interpolation on a matrix.
    // Construct with a vector of x1 values, a vector of x2 values, and a
    // matrix of tabulated function values yij. Then call interp for
    // interpolated values.
    protected int m, n;
    protected final double[][] y;
    protected final double[] x1;
    protected final double[] yv;
    protected List<Spline_interp> srp;

    public Spline2D_interp(final double[] x1v, final double[] x2v, final double[][] ym) throws NRException {
        m = (x1v.length);
        n = (x2v.length);
        y = doub_mat(ym);
        yv = doub_arr(m);
        x1 = (x1v);
        srp = new ArrayList<Spline_interp>();
        for (int i = 0; i < m; i++)
            srp.add(new Spline_interp(x2v, y[i]));
        // Save an array of pointers to 1-dim row splines.
    }

    // ~Spline2D_interp(){
    // for (Int i=0;i<m;i++) delete srp[i]; We need a destructor to clean up.
    // }

    public double interp(final double x1p, final double x2p) throws NRException {
        for (int i = 0; i < m; i++)
            yv[i] = srp.get(i).interp(x2p);
        // Interpolate on each row.
        Spline_interp scol = 
            new Spline_interp(x1, yv); // Construct the column
                                          // spline,
        return scol.interp(x1p); // and evaluate it.
    }

}
