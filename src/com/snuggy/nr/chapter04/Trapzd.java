
package com.snuggy.nr.chapter04;

import com.snuggy.nr.util.*;

public class Trapzd<T extends Func_Doub_To_Doub> extends Quadrature {

    // Routine implementing the extended trapezoidal rule.

    private double a, b, s; // Limits of integration and current value of
                            // integral.
    private T func;

    public Trapzd() {
    }

    public Trapzd(final T funcc, final double aa, final double bb) {
        func = (funcc);
        a = (aa);
        b = (bb);
        n = 0;
    }

    // The constructor takes as inputs func, the function or functor to be
    // integrated between limits a and b, also input.

    public double next() throws NRException {
        // Returns the nth stage of refinement of the extended trapezoidal
        // rule. On the first call (n=1), the routine returns the crudest
        // estimate of
        // R b
        // a f.x/dx. Subsequent calls set n=2,3,... and improve the accuracy
        // by adding 2n-2 additional interior points.
        double x, tnm, sum, del;
        int it, j;
        n++;
        if (n == 1) {
            return (s = 0.5 * (b - a) * (func.eval(a) + func.eval(b)));
        } else {
            for (it = 1, j = 1; j < n - 1; j++)
                it <<= 1;
            tnm = it;
            del = (b - a) / tnm; // This is the spacing of the points to be
                                 // added.
            x = a + 0.5 * del;
            for (sum = 0.0, j = 0; j < it; j++, x += del)
                sum += func.eval(x);
            s = 0.5 * (s + (b - a) * sum / tnm); // This replaces s by its
                                                 // refined value.
            return s;
        }
    }

}
