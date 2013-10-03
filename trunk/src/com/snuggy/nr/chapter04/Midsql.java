
package com.snuggy.nr.chapter04;

import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Midsql<T extends Func_Doub_To_Doub> extends Midpnt<T> {

    // This routine is an exact replacement for midpnt, except that
    // it allows for an inverse square-root singularity in the integrand
    // at the lower limit aa.

    private double aorig;

    public double func(final double x) throws NRException {
        return 2.0 * x * super.funk.eval(aorig + x * x); // Effect the change of
                                                         // variable.
    }

    public Midsql(final T funcc, final double aa, final double bb) {
        super(funcc, aa, bb);
        aorig = (aa);
        super.a = 0;
        super.b = sqrt(bb - aa);
    }

}
