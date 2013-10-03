
package com.snuggy.nr.chapter03;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.chapter02.*;
import com.snuggy.nr.util.*;

public class Krig<T extends Func_Doub_To_Doub> {

    // Object for interpolation by kriging, using npt points in ndim
    // dimensions. Call constructor once, then interp as many times as desired.
    private final double[][] x;
    private final T vgram;
    private int ndim, npt;
    private double lastval, lasterr; // Most recently computed value and (if com
    double[] y, dstar, vstar, yvi; // puted) error.
    private double[][] v;
    private LUdcmp vi;

    public Krig(final double[][] xx, final double[] yy, final T vargram) throws NRException {
        this(xx, yy, vargram, null, 0);
    }

    public Krig(final double[][] xx, final double[] yy, final T vargram, 
                final double[] err_arr, final int err_off)
            throws NRException {
        // Constructor. The npt  ndim matrix xx inputs the data points, the
        // vector yy the function values. vargram is the variogram function
        // or functor. The argument err is not used for interpolation; see
        // 15.9.

        x = (xx);
        vgram = (vargram);
        npt = (nrows(xx));
        ndim = (ncols(xx));
        dstar = doub_arr(npt + 1);
        vstar = doub_arr(npt + 1);
        v = doub_mat(npt + 1, npt + 1);
        y = doub_arr(npt + 1);
        yvi = doub_arr(npt + 1);
        int i, j;
        for (i = 0; i < npt; i++) { // Fill Y and V.
            y[i] = yy[i];
            for (j = i; j < npt; j++) {
                v[i][j] = v[j][i] = vgram.eval(rdist(x[i], 0, x[j], 0));
            }
            v[i][npt] = v[npt][i] = 1.;
        }
        v[npt][npt] = y[npt] = 0.;
        if (err_arr != null)
            for (i = 0; i < npt; i++)
                v[i][i] -= SQR(err_arr[err_off + i]); // 15.9.
        vi = new LUdcmp(v);
        vi.solve(y, yvi);
    }

    // ~Krig() { delete vi; }

    public double interp(final double[] xstar) throws NRException {
        // Return an interpolated value at the point xstar.
        int i;
        for (i = 0; i < npt; i++)
            vstar[i] = vgram.eval(rdist(xstar, 0, x[i], 0));
        vstar[npt] = 1.;
        lastval = 0.;
        for (i = 0; i <= npt; i++)
            lastval += yvi[i] * vstar[i];
        return lastval;
    }

    public double interp(final double[] xstar, final double esterr_ref[]) 
            throws NRException {
        // Return an interpolated value at the point xstar, and return its
        // estimated error as esterr.
        lastval = interp(xstar);
        vi.solve(vstar, dstar);
        lasterr = 0;
        for (int i = 0; i <= npt; i++)
            lasterr += dstar[i] * vstar[i];
        esterr_ref[0] = lasterr = sqrt(MAX(0., lasterr));
        return lastval;
    }

    private double rdist(final double[] x1_arr, final int x1_off, 
                            final double[] x2_arr, final int x2_off) {
        // Utility used internally. Cartesian distance between two points.
        double d = 0.0;
        for (int i = 0; i < ndim; i++)
            d += SQR(x1_arr[x1_off + i] - x2_arr[x2_off + i]);
        return sqrt(d);
    }

}
