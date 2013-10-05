
package com.snuggy.nr.chapter03;

import static com.snuggy.nr.refs.Refs.*;
import static com.snuggy.nr.util.Static.*;

import com.snuggy.nr.util.*;

public class BaryRat_interp extends Base_interp {

    // Barycentric rational interpolation object. After constructing the
    // object, call interp for interpolated values. Note that no error
    // estimate dy is calculated.

    private final double[] w;
    private int d;

    // BaryRat_interp(VecDoub_I &xv, VecDoub_I &yv, Int dd);
    // Doub rawinterp(Int jl, Doub x);
    // Doub interp(Doub x);

    public BaryRat_interp(final double[] xv, final double[] yv, final int dd) throws NRException {
        // Constructor arguments are x and y vectors of length n, and order d of
        // desired approximation.
        super(xv, $(yv, 0), xv.length);
        w = doub_arr(n);
        d = (dd);

        if (n <= d)
            throw new NRException("d too large for number of points in BaryRat_interp");
        for (int k = 0; k < n; k++) { // Compute weights from equation (3.4.10).
            int imin = MAX(k - d, 0);
            int imax = k >= n - d ? n - d - 1 : k;
            double temp = ((imin & 1) != 0) ? -1.0 : 1.0;
            double sum = 0.0;
            for (int i = imin; i <= imax; i++) {
                int jmax = MIN(i + d, n - 1);
                double term = 1.0;
                for (int j = i; j <= jmax; j++) {
                    if (j == k)
                        continue;
                    term *= (xx.$(k) - xx.$(j));
                }
                term = temp / term;
                temp = -temp;
                sum += term;
            }
            w[k] = sum;
        }
    }

    public double rawinterp(final int jl, final double x) throws NRException {
        // Use equation (3.4.9) to compute the barycentric rational
        // interpolant. Note that jl is not used since the approximation
        // is global; it is included only for compatibility with Base_interp.

        double num = 0, den = 0;
        for (int i = 0; i < n; i++) {
            double h = x - xx.$(i);
            if (h == 0.0) {
                return yy.$(i);
            } else {
                double temp = w[i] / h;
                num += temp * yy.$(i);
                den += temp;
            }
        }
        return num / den;
    }

    public double interp(final double x) throws NRException {
        // No need to invoke hunt or locate since the interpolation is
        // global, so override interp to simply call rawinterp directly
        // with a dummy value of jl.
        return rawinterp(1, x);
    }

}
