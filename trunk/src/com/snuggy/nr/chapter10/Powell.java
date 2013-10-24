
package com.snuggy.nr.chapter10;

import static java.lang.Math.*;
import static com.snuggy.nr.util.Static.*;

import com.snuggy.nr.util.*;

import static com.snuggy.nr.refs.Refs.*;

public class Powell<T extends Func_DoubVec_To_Doub> extends Linemethod<T> {

    // Multidimensional minimization by Powell’s method.
    private int iter;
    protected double fret; // Value of the function at the minimum.

    // using Linemethod<T>::func; Variables from a templated base class are not
    // autousing
    // Linemethod<T>::linmin; matically inherited.
    // using Linemethod<T>::p;
    // using Linemethod<T>::xi;

    private final double ftol;

    public Powell(final T func) throws NRException {
        this(func, 3.0e-8);
    }

    public Powell(final T func, final double ftoll) throws NRException {
        // Constructor arguments are func, the function or functor to be
        // minimized, and an optional argument ftoll, the fractional tolerance
        // in the function value such that failure to decrease by more than this
        // amount on one iteration signals doneness.
        super(func);
        ftol = (ftoll);
    }

    public double[] minimize(final double[] pp) throws NRException {
        // Minimization of a function or functor n variables. Input consists of
        // an initial starting point pp[0..n-1]. The initial matrix
        // ximat[0..n-1][0..n-1], whose columns contain the initial set of
        // directions, is set to the n unit vectors. Returned is the best point
        // found, at which point fret is the minimum function value and iter is
        // the number of iterations taken.
        int n = pp.length;
        double[][] ximat = doub_mat(n, n, 0.0);
        for (int i = 0; i < n; i++)
            ximat[i][i] = 1.0;
        return minimize(pp, ximat);
    }

    public double[] minimize(final double[] pp, final double[][] ximat) throws NRException {
        // Alternative interface: Input consists of the initial starting point
        // pp[0..n-1] and an initial matrix ximat[0..n-1][0..n-1], whose columns
        // contain the initial set of directions. On output ximat is the
        // then-current direction set.
        final int ITMAX = 200; // Maximum allowed iterations.
        final double TINY = 1.0e-25; // A small number.
        double fptt;
        int n = pp.length;
        $$(p, pp);
        double[] pt = doub_vec(n), ptt = doub_vec(n);
        //xi.resize(n);
        $$(xi, doub_vec(n));
        fret = func.eval(p.$());
        for (int j = 0; j < n; j++)
            pt[j] = p.$()[j]; // Save the initial point.
        for (iter = 0;; ++iter) {
            double fp = fret;
            int ibig = 0;
            double del = 0.0; // Will be the biggest function decrease.
            for (int i = 0; i < n; i++) { // In each iteration, loop over all
                                          // directions in the set.
                for (int j = 0; j < n; j++)
                    xi.$()[j] = ximat[j][i]; // Copy the direction,
                fptt = fret;
                fret = linmin(); // minimize along it,
                if (fptt - fret > del) { // and record it if it is the largest
                                         // decrease
                    del = fptt - fret; // so far.
                    ibig = i + 1;
                }
            } // Here comes the termination criterion:
            if (2.0 * (fp - fret) <= ftol * (abs(fp) + abs(fret)) + TINY) {
                return p.$$();
            }
            if (iter == ITMAX)
                throw new NRException("powell exceeding maximum iterations.");
            for (int j = 0; j < n; j++) { // Construct the extrapolated point
                                          // and the
                // average direction moved. Save the old starting point.
                ptt[j] = 2.0 * p()[j] - pt[j];
                xi()[j] = p()[j] - pt[j];
                pt[j] = p()[j];
            }
            fptt = func.eval(ptt); // Function value at extrapolated point.
            if (fptt < fp) {
                double t = 2.0 * (fp - 2.0 * fret + fptt) * SQR(fp - fret - del) - del * SQR(fp - fptt);
                if (t < 0.0) {
                    fret = linmin(); // Move to the minimum of the new direc
                    for (int j = 0; j < n; j++) { // tion, and save the new
                                                  // direction.
                        ximat[j][ibig - 1] = ximat[j][n - 1];
                        ximat[j][n - 1] = xi()[j];
                    }
                }
            }
        }
    }

}
