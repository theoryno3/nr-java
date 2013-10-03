
package com.snuggy.nr.chapter03;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Powvargram implements Func_Doub_To_Doub {

    // Functor for variogram v.r/ D ?r?, where ? is specified, ? is fitted
    // from the data.

    private double alph, bet, nugsq;

    public Powvargram(final double[][] x, final double[] y) {
        this(x, y, 1.5, 0.0);
    }

    public Powvargram(final double[][] x, final double[] y, final double beta) {
        this(x, y, beta, 0.0);
    }

    public Powvargram(final double[][] x, final double[] y, final double beta, final double nug) {
        // Constructor. The npt  ndim matrix x inputs the data points, the
        // vector y the function values, beta the value of ?. For interpolation,
        // the default value of beta is usually adequate. For the (rare) use
        // of nug see 15.9.

        bet = (beta);
        nugsq = (nug * nug);
        int i, j, k, npt = nrows(x), ndim = ncols(x);
        double rb, num = 0.0, denom = 0.0;
        for (i = 0; i < npt; i++)
            for (j = i + 1; j < npt; j++) {
                rb = 0.0;
                for (k = 0; k < ndim; k++)
                    rb += SQR(x[i][k] - x[j][k]);
                rb = pow(rb, 0.5 * beta);
                num += rb * (0.5 * SQR(y[i] - y[j]) - nugsq);
                denom += SQR(rb);
            }
        alph = num / denom;
    }

    public double eval(final double r) {
        return nugsq + alph * pow(r, bet);
    }

}
