package com.snuggy.nr.chapter09;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class NRfdjac<T extends Func_DoubVec_To_DoubVec> {

    // Computes forward-difference approximation to Jacobian.
    private final double EPS; // Set to approximate square root of the machine
                              // pre
    private final T func; // cision.

    public NRfdjac(final T funcc) {
        // Initialize with user-supplied function or functor that returns the
        // vector of functions to be zeroed.
        EPS = (1.0e-8);
        func = (funcc);
    }

    public double[][] eval(final double[] x, final double[] fvec) {
        // Returns the Jacobian array df[0..n-1][0..n-1]. On input, x[0..n-1]
        // is the point at which the Jacobian is to be evaluated and
        // fvec[0..n-1] is the vector of function values at the point.
        final int n = x.length;
        final double[][] df = doub_mat(n, n);
        final double[] xh = doub_vec(x);
        for (int j = 0; j < n; j++) {
            final double temp = xh[j];
            double h = EPS * abs(temp);
            if (h == 0.0)
                h = EPS;
            xh[j] = temp + h; // Trick to reduce finite-precision error.
            h = xh[j] - temp;
            double[] f = func.eval(xh);
            xh[j] = temp;
            for (int i = 0; i < n; i++)
                // Forward difference formula.
                df[i][j] = (f[i] - fvec[i]) / h;
        }
        return df;
    }

}
