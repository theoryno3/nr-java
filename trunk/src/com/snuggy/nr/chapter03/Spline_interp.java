
package com.snuggy.nr.chapter03;

import static com.snuggy.nr.util.Static.*;

import com.snuggy.nr.util.*;

public class Spline_interp extends Base_interp {

    // Cubic spline interpolation object. Construct with x and y vectors,
    // and (optionally) values of the first derivative at the endpoints, then
    // call interp for interpolated values.

    private double[] y2;

    public Spline_interp() {
    }
    
    public Spline_interp(final double[] xv, final double[] yv_arr, final int yv_off) {
        this(xv, yv_arr, yv_off, 1.e99, 1.e99);
    }

    public Spline_interp(final double[] xv, final double[] yv_arr, final int yv_off, final double yp1) {
        this(xv, yv_arr, yv_off, yp1, 1.e99);
    }

    public Spline_interp(final double[] xv, final double[] yv_arr, final int yv_off, final double yp1, final double ypn) {
        super(xv, yv_arr, yv_off, 2);
        y2 = doub_arr(xv.length);
        sety2(xv, 0, yv_arr, yv_off, yp1, ypn);
    }

    // public Spline_interp(double[] xv, final double yv_arr, final int yv_off,
    // double yp1=1.e99, Doub ypn=1.e99)
    // : Base_interp(xv,yv,2), y2(xv.size())
    // {sety2(&xv[0],yv,yp1,ypn);}

    // public void sety2(final double xv_arr, final int xv_off, final double
    // yv_arr, final int yv_off, double yp1, double ypn);
    // Doub rawinterp(Int jl, Doub xv);

    // For now, you can ignore the second constructor; it will be used later
    // for two-dimensional spline interpolation. The user interface differs
    // from previous ones only in the addition of two constructor arguments,
    // used to set the values of the first derivatives at the endpoints, y00
    // and y0 N1. These are coded with default values that signal that you
    // want a natural spline, so they can be omitted in most situations. Both
    // constructors invoke sety2 to do the actual work of computing, and
    // storing, the second derivatives.

    public void sety2(final double[] xv_arr, final int xv_off, 
            final double[] yv_arr, final int yv_off, 
            final double yp1, final double ypn) {
        // This routine stores an array y2[0..n-1] with second derivatives of
        // the interpolating function at the tabulated points pointed to by xv,
        // using function values pointed to by yv. If yp1 and/or ypn are equal
        // to 1  1099 or larger, the routine is signaled to set the
        // corresponding boundary condition for a natural spline, with zero
        // second derivative on that boundary; otherwise, they are the values
        // of the first derivatives at the endpoints.

        int i, k;
        double p, qn, sig, un;
        int n = y2.length;
        double[] u = doub_arr(n - 1);
        if (yp1 > 0.99e99) // The lower boundary condition is set either to be “
            y2[0] = u[0] = 0.0; // natural”
        else { // or else to have a specified first derivative.
            y2[0] = -0.5;
            u[0] = (3.0 / (xv_arr[xv_off + 1] - xv_arr[xv_off + 0]))
                    * ((yv_arr[yv_off + 1] - yv_arr[yv_off + 0]) / 
                            (xv_arr[xv_off + 1] - xv_arr[xv_off + 0]) - yp1);
        }
        for (i = 1; i < n - 1; i++) { // This is the decomposition loop of the
                                      // tridiagonal algorithm.
            // y2 and u are used for temporary
            // storage of the decomposed
            // factors.
            sig = (xv_arr[xv_off + i] - xv_arr[xv_off + i - 1]) / 
                    (xv_arr[xv_off + i + 1] - xv_arr[xv_off + i - 1]);
            p = sig * y2[i - 1] + 2.0;
            y2[i] = (sig - 1.0) / p;
            u[i] = (yv_arr[yv_off + i + 1] - yv_arr[yv_off + i]) / 
                    (xv_arr[xv_off + i + 1] - xv_arr[xv_off + i])
                    - (yv_arr[yv_off + i] - yv_arr[yv_off + i - 1]) / 
                    (xv_arr[xv_off + i] - xv_arr[xv_off + i - 1]);
            u[i] = (6.0 * u[i] / 
                    (xv_arr[xv_off + i + 1] - xv_arr[xv_off + i - 1]) - 
                    sig * u[i - 1]) / p;
        }
        if (ypn > 0.99e99) // The upper boundary condition is set either to be
            qn = un = 0.0; // “natural”
        else { // or else to have a specified first derivative.
            qn = 0.5;
            un = (3.0 / (xv_arr[xv_off + n - 1] - xv_arr[xv_off + n - 2]))
                    * (ypn - (yv_arr[yv_off + n - 1] - yv_arr[yv_off + n - 2])
                        / (xv_arr[xv_off + n - 1] - xv_arr[xv_off + n - 2]));
        }
        y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);
        for (k = n - 2; k >= 0; k--)
            // This is the backsubstitution loop of the tridiagonal
            y2[k] = y2[k] * y2[k + 1] + u[k]; // algorithm.
    }

    // Note that, unlike the previous object Poly_interp, Spline_interp stores
    // data that depend on the contents of your array of yi ’s at its time of
    // creation — a whole vector y2. Although we didn’t point it out, the
    // previous interpolation object actually allowed the misuse of altering
    // the contents of their x and y arrays on the fly (as long as the lengths
    // didn’t change). If you do that with Spline_interp, you’ll get definitely
    // wrong answers! The required rawinterp method, never called directly by
    // the users, uses the stored y2 and implements equation (3.3.3):

    public double rawinterp(final int jl, final double x) throws NRException {
        // Given a value x, and using pointers to data xx and yy, and the
        // stored vector of second derivatives y2, this routine returns the
        // cubic spline interpolated value y.

        int klo = jl, khi = jl + 1;
        double y, h, b, a;
        h = xx_arr[xx_off + khi] - xx_arr[xx_off + klo];
        if (h == 0.0)
            throw new NRException("Bad input to routine splint"); // The xa’s
                                                                  // must be dis
        a = (xx_arr[xx_off + khi] - x) / h; // tinct.
        b = (x - xx_arr[xx_off + klo]) / h; // Cubic spline polynomial is now
                                            // evaluated.
        y = a * yy_arr[yy_off + klo] + b * yy_arr[yy_off + khi] + 
                ((a * a * a - a) * y2[klo] + (b * b * b - b) * y2[khi]) * 
                (h * h) / 6.0;
        return y;
    }

}
