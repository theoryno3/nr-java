package com.snuggy.nr.chapter10;

import static java.lang.Math.*;

import com.snuggy.nr.util.*;

import static com.snuggy.nr.refs.Refs.*;

import com.snuggy.nr.refs.*;

public class Golden extends Bracketmethod {

    // Golden section search for minimum.
    private double xmin;
    @SuppressWarnings("unused")
    private double fmin;
    private final double tol;

    public Golden() {
        this(3.0e-8);
    }

    public Golden(final double toll) {
        tol = (toll);
    }

    public <T extends Func_Doub_To_Doub> double minimize(final T func) throws NRException {
        // Given a function or functor f, and given a bracketing triplet of
        // abscissas ax, bx, cx (such that bx is between ax and cx, and f(bx)
        // is less than both f(ax) and f(cx)), this routine performs a golden
        // section search for the minimum, isolating it to a fractional
        // precision of about tol. The abscissa of the minimum is returned as
        // xmin, and the function value at the minimum is returned as min, the
        // returned function value.
        final double R = 0.61803399, C = 1.0 - R; // The golden ratios.
        $double x1 = $(0.0), x2 = $(0.0);
        $double x0 = $(ax.$()); // At any given time we will keep track of four
        $double x3 = $(cx.$()); // points, x0,x1,x2,x3.
        if (abs(cx.$() - bx.$()) > abs(bx.$() - ax.$())) { // Make x0 to x1 the
                                                           // smaller segment,
            $(x1, bx);
            x2.$(bx.$() + C * (cx.$() - bx.$())); // and fill in the new point
                                                  // to be tried.
        } else {
            $(x2, bx);
            x1.$(bx.$() - C * (bx.$() - ax.$()));
        }
        $double f1 = $(func.eval(x1.$())); // The initial function evaluations.
                                           // Note that
        // we never need to evaluate the function at the original endpoints.
        $double f2 = $(func.eval(x2.$()));
        while (abs(x3.$() - x0.$()) > tol * (abs(x1.$()) + abs(x2.$()))) {
            if (f2.$() < f1.$()) { // One possible outcome,
                shft3(x0, x1, x2, R * x2.$() + C * x3.$()); // its housekeeping,
                shft2(f1, f2, func.eval(x2.$())); // and a new function
                                                  // evaluation.
            } else { // The other outcome,
                shft3(x3, x2, x1, R * x1.$() + C * x0.$());
                shft2(f2, f1, func.eval(x1.$())); // and its new function
                                                  // evaluation.
            }
        } // Back to see if we are done.
        if (f1.$() < f2.$()) { // We are done. Output the best of the two
            xmin = x1.$(); // current values.
            fmin = f1.$();
        } else {
            xmin = x2.$();
            fmin = f2.$();
        }
        return xmin;
    }

}
