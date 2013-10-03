
package com.snuggy.nr.chapter04;

import com.snuggy.nr.util.*;

public class Midinf<T extends Func_Doub_To_Doub> extends Midpnt<T> {

    // This routine is an exact replacement for midpnt, i.e., returns the
    // nth stage of refinement of the integral of funcc from aa to bb, except
    // that the function is evaluated at evenly spaced points in 1=x rather
    // than in x. This allows the upper limit bb to be as large and positive
    // as the computer allows, or the lower limit aa to be as large and
    // negative, but not both. aa and bb must have the same sign.

    public double func(final double x) throws NRException {
        return super.funk.eval(1.0 / x) / (x * x); // Effect the change of
                                                   // variable.
    }

    public Midinf(final T funcc, final double aa, final double bb) {
        super(funcc, aa, bb);
        super.a = 1.0 / bb; // Set the limits of integration.
        super.b = 1.0 / aa;
    }

}
