package com.snuggy.nr.chapter17;

import static com.snuggy.nr.util.Static.*;

import com.snuggy.nr.chapter02.*;
import com.snuggy.nr.util.*;

public class Static {

    public static void rk4(final double[] y, final double[] dydx, final double x, final double h, final double[] yout,
            Func_Doub_DoubArr_DoubArr_To_Void derivs) {
        // Given values for the variables y[0..n-1] and their derivatives
        // dydx[0..n-1] known at x, use the fourth-order Runge-Kutta method to
        // advance the solution over an interval h and return the incremented
        // variables as yout[0..n-1]. The user supplies the routine
        // derivs(x,y,dydx), which returns derivatives dydx at x.
        int n = y.length;
        final double[] dym = doub_arr(n), dyt = doub_arr(n), yt = doub_arr(n);
        double hh = h * 0.5;
        double h6 = h / 6.0;
        double xh = x + hh;
        for (int i = 0; i < n; i++)
            yt[i] = y[i] + hh * dydx[i]; // First step.
        derivs.eval(xh, yt, dyt); // Second step.
        for (int i = 0; i < n; i++)
            yt[i] = y[i] + hh * dyt[i];
        derivs.eval(xh, yt, dym); // Third step.
        for (int i = 0; i < n; i++) {
            yt[i] = y[i] + h * dym[i];
            dym[i] += dyt[i];
        }
        derivs.eval(x + h, yt, dyt); // Fourth step.
        for (int i = 0; i < n; i++)
            // Accumulate increments with
            yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]); // proper
                                                                     // weights.
    }

    public static void sparmatfill(NRsparseCol[] sparmat, final double[][] fullmat) throws NRException {
        // Utility that fills a sparse matrix from a full one. See 2.7.
        int n, m, nz, nn = nrows(fullmat), mm = ncols(fullmat);
        if (sparmat.length != mm)
            throw new NRException("bad sizes");
        for (m = 0; m < mm; m++) {
            for (nz = n = 0; n < nn; n++)
                if (fullmat[n][m] != 0.0)
                    nz++;
            sparmat[m].resize(nn, nz);
            for (nz = n = 0; n < nn; n++)
                if (fullmat[n][m] != 0.0) {
                    sparmat[m].row_ind()[nz] = n;
                    sparmat[m].val()[nz++] = fullmat[n][m];
                }
        }
    }
}
