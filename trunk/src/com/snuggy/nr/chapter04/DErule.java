
package com.snuggy.nr.chapter04;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class DErule<T extends Func_Doub_Doub_To_Doub> extends Quadrature {

    // Structure for implementing the DE rule.
    private double a, b, hmax, s;
    private T func;

    public DErule(final T funcc, final double aa, final double bb) {
        this(funcc, aa, bb, 3.7);
    }

    public DErule(final T funcc, final double aa, final double bb, final double hmaxx) {
        func = (funcc);
        a = (aa);
        b = (bb);
        hmax = (hmaxx);
        n = 0;
    }

    // Constructor. funcc is the function or functor that provides the
    // function to be integrated between limits aa and bb, also input. The
    // function operator in funcc takes two arguments, x and i, as described
    // in the text. The range of integration in the transformed variable t
    // is .hmaxx; hmaxx/. Typical values of hmaxx are 3.7 for logarithmic or
    // milder singularities, and 4.3 for square-root singularities, as
    // discussed in the text.

    @Override
    public double next() {
        // On the first call to the function next (n D 1), the routine returns
        // the crudest estimate of R b a f.x/dx. Subsequent calls to next
        // (n D 2;3; : : :) will improve the accuracy by adding 2n1 additional
        // interior points.
        double del, fact, q, sum, t, twoh;
        int it, j;
        n++;
        if (n == 1) {
            fact = 0.25;
            return s = hmax * 2.0 * (b - a) * fact * 
                        func.eval(0.5 * (b + a), 0.5 * (b - a));
        } else {
            for (it = 1, j = 1; j < n - 1; j++)
                it <<= 1;
            twoh = hmax / it; // Twice the spacing of the points to be added.
            t = 0.5 * twoh;
            for (sum = 0.0, j = 0; j < it; j++) {
                q = exp(-2.0 * sinh(t));
                del = (b - a) * q / (1.0 + q);
                fact = q / SQR(1.0 + q) * cosh(t);
                sum += fact * (func.eval(a + del, del) + func.eval(b - del, del));
                t += twoh;
            }
            return s = 0.5 * s + (b - a) * twoh * sum; // Replace s by its
                                                       // refined value and
                                                       // return.
        }
    }

}
