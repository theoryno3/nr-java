
package com.snuggy.nr.chapter10;

import static com.snuggy.nr.refs.Refs.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.refs.*;
import com.snuggy.nr.util.*;

public class Bracketmethod {

    // Base class for one-dimensional minimization routines. Provides a
    // routine to bracket a minimum and several utility functions.
    protected $double ax = $(0.0);
    protected $double bx = $(0.0);
    protected $double cx = $(0.0);
    protected $double fa = $(0.0);
    protected $double fb = $(0.0);
    protected $double fc = $(0.0);
    
    public double ax() {
        return ax.$();
    }
    
    public double bx() {
        return bx.$();
    }
    
    public double cx() {
        return cx.$();
    }
    
    public double fa() {
        return fa.$();
    }

    public double fb() {
        return fb.$();
    }

    public double fc() {
        return fc.$();
    }

    public <T extends Func_Doub_To_Doub> void bracket(final double a, final double b, final T func) throws NRException {
        // Given a function or functor func, and given distinct initial points
        // ax and bx, this routine searches in the downhill direction (defined
        // by the function as evaluated at the initial points) and returns new
        // points ax, bx, cx that bracket a minimum of the function. Also
        // returned are the function values at the three points, fa, fb, and fc.
        final double GOLD = 1.618034, GLIMIT = 100.0, TINY = 1.0e-20;
        // Here GOLD is the default ratio by which successive intervals are
        // magnified and GLIMIT is the maximum magnification allowed for a
        // parabolic-fit step.
        $(ax, a);
        $(bx, b);
        $double fu = $(0.0);
        $(fa, func.eval(ax.$()));
        $(fb, func.eval(bx.$()));
        if (fb.$() > fa.$()) { // Switch roles of a and b so that we can
                                     // go
            // SWAP(ax,bx); // downhill in the direction from a to b.
            SWAP(ax, bx); // downhill in the direction from a to b.
            // SWAP(fb,fa);
            SWAP(fb, fa);
        }
        $(cx, bx.$() + GOLD * (bx.$() - ax.$())); // First guess
                                                                // for c.
        $(fc, func.eval(cx.$()));
        while (fb.$() > fc.$()) { // Keep returning here until we bracket.
            double r = (bx.$() - ax.$()) * (fb.$() - fc.$()); // Compute
                                                                          // u
                                                                          // by
                                                                          // parabolic
                                                                          // extrapolation
                                                                          // from
            // a; b; c. TINY is used to prevent any possible division by zero.
            double q = (bx.$() - cx.$()) * (fb.$() - fa.$());
            $double u = $(bx.$() - ((bx.$() - cx.$()) * q - (bx.$() - ax.$()) * r)
                    / (2.0 * SIGN(MAX(abs(q - r), TINY), q - r)));
            double ulim = bx.$() + GLIMIT * (cx.$() - bx.$());
            // We won’t go farther than this. Test various possibilities:
            if ((bx.$() - u.$()) * (u.$() - cx.$()) > 0.0) { // Parabolic
                                                                         // u is
                                                                         // between
                                                                         // b
                                                                         // and
                                                                         // c:
                                                                         // try
                                                                         // it.
                $(fu, func.eval(u.$()));
                if (fu.$() < fc.$()) { // Got a minimum between b and c.
                    $(ax, bx);
                    $(bx, u);
                    $(fa, fb);
                    $(fb, fu);
                    return;
                } else if (fu.$() > fb.$()) { // Got a minimum between
                                                    // between a and u.
                    $(cx, u);
                    $(fc, fu);
                    return;
                }
                u.$(cx.$() + GOLD * (cx.$() - bx.$())); // Parabolic
                                                                       // fit
                                                                       // was no
                                                                       // use.
                                                                       // Use
                                                                       // default
                                                                       // mag
                $(fu, func.eval(u.$())); // nification.
            } else if ((cx.$() - u.$()) * (u.$() - ulim) > 0.0) { // Parabolic
                                                                           // fit
                                                                           // is
                                                                           // between
                                                                           // c
                                                                           // and
                $(fu, func.eval(u.$())); // its allowed limit.
                if (fu.$() < fc.$()) {
                    shft3(bx, cx, u, u.$() + GOLD * (u.$() - cx.$()));
                    shft3(fb, fc, fu, func.eval(u.$()));
                }
            } else if ((u.$() - ulim) * (ulim - cx.$()) >= 0.0) { // Limit
                                                                        // parabolic
                                                                        // u to
                                                                        // maximum
                u.$(ulim); // allowed value.
                $(fu, func.eval(u.$()));
            } else { // Reject parabolic u, use default magnifica
                $(u, cx.$() + GOLD * (cx.$() - bx.$())); // tion.
                $(fu, func.eval(u.$()));
            }
            shft3(ax, bx, cx, u.$()); // Eliminate oldest point
                                                     // and continue.
            shft3(fa, fb, fc, fu.$());
        }
    }

    private void SWAP(final $double x, final $double y) {
        double t = x.$();
        $(x, y);
        $(y, t);
    }

    public void shft2(final $double a, final $double b, final double c) {
        // Utility function used in this structure or others derived from it.
        a.$(b.$());
        b.$(c);
    }

    public void shft3(final $double a, final $double b, final $double c, final double d) {
        $(a, b);
        $(b, c);
        $(c, d);
    }

    public void mov3(final $double a, final $double b, final $double c, final double d, final double e,
            final double f) {
        a.$(d);
        b.$(e);
        c.$(f);
    }

}
