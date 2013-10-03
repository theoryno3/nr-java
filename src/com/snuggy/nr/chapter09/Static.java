
package com.snuggy.nr.chapter09;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Static {

    public static <T extends Func_Doub_To_Doub> double zbrent(final T func, final double x1, final double x2, final double tol)
            throws NRException {
        // Using Brent’s method, return the root of a function or functor func
        // known to lie between x1 and x2. The root will be refined until its
        // accuracy is tol.
        final int ITMAX = 100; // Maximum allowed number of iterations.
        final double EPS = EPS(); // numeric_limits<Doub>::epsilon();
        // Machine floating-point precision.
        double a = x1, b = x2, c = x2, d = 0.0, e = 0.0, fa = func.eval(a), fb = func.eval(b), fc, p, q, r, s, tol1, xm;
        if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
            throw new NRException("Root must be bracketed in zbrent");
        fc = fb;
        for (int iter = 0; iter < ITMAX; iter++) {
            if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
                c = a; // Rename a, b, c and adjust bounding interval
                fc = fa; // d.
                e = d = b - a;
            }
            if (abs(fc) < abs(fb)) {
                a = b;
                b = c;
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }
            tol1 = 2.0 * EPS * abs(b) + 0.5 * tol; // Convergence check.
            xm = 0.5 * (c - b);
            if (abs(xm) <= tol1 || fb == 0.0)
                return b;
            if (abs(e) >= tol1 && abs(fa) > abs(fb)) {
                s = fb / fa; // Attempt inverse quadratic interpolation.
                if (a == c) {
                    p = 2.0 * xm * s;
                    q = 1.0 - s;
                } else {
                    q = fa / fc;
                    r = fb / fc;
                    p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                    q = (q - 1.0) * (r - 1.0) * (s - 1.0);
                }
                if (p > 0.0)
                    q = -q; // Check whether in bounds.
                p = abs(p);
                double min1 = 3.0 * xm * q - abs(tol1 * q);
                double min2 = abs(e * q);
                if (2.0 * p < (min1 < min2 ? min1 : min2)) {
                    e = d; // Accept interpolation.
                    d = p / q;
                } else {
                    d = xm; // Interpolation failed, use bisection.
                    e = d;
                }
            } else { // Bounds decreasing too slowly, use bisection.
                d = xm;
                e = d;
            }
            a = b; // Move last best guess to a.
            fa = fb;
            if (abs(d) > tol1) // Evaluate new trial root.
                b += d;
            else
                b += SIGN(tol1, xm);
            fb = func.eval(b);
        }
        throw new NRException("Maximum number of iterations exceeded in zbrent");
    }

}
