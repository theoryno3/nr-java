package com.snuggy.nr.chapter03;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public abstract class Base_interp {

    // Abstract base class used by all interpolation routines in this chapter.
    // Only the routine interp is called directly by the user.

    protected final int n, mm, dj;
    protected final double[] xx_arr;
    protected final int xx_off;
    protected double[] yy_arr;
    protected int yy_off;
    protected int cor, jsav;
    
    public Base_interp() {
	    n = 0; mm = 0; dj = 0;
	    xx_arr = null;
	    xx_off = 0;
	    yy_arr = null;
	    yy_off = 0;
    }

    public Base_interp(final double[] x, final double[] y_arr, final int y_off, final int m) {
        // Constructor: Set up for interpolating on a table of x’s and y’s
        // of length m. Normally called by a derived class, not by the user.
        n = (x.length);
        mm = (m);
        jsav = (0);
        cor = (0);
        // xx = (&x[0]);
        xx_arr = x;
        xx_off = 0;
        yy_arr = (y_arr);
        yy_off = y_off;
        dj = MIN(1, (int) pow((double) n, 0.25));
    }

    public double interp(final double x) throws NRException {
        // Given a value x, return an interpolated value, using data
        // pointed to by xx and yy.
        int jlo = (cor != 0) ? hunt(x) : locate(x);
        return rawinterp(jlo, x);
    }

    // Int locate(final double x); See definitions below.
    // Int hunt(final double x);

    public abstract double rawinterp(final int jlo, final double x) throws NRException;

    // Derived classes provide this as the actual interpolation method.

    public int locate(final double x) throws NRException {
        // Given a value x, return a value j such that x is (insofar as
        // possible) centered in the subrange xx[j..j+mm-1], where xx is
        // the stored pointer. The values in xx must be monotonic, either
        // increasing or decreasing. The returned value is not less than 0,
        // nor greater than n-1.

        int ju, jm, jl;
        if (n < 2 || mm < 2 || mm > n)
            throw new NRException("locate size error");
        boolean ascnd = (xx_arr[xx_off + n - 1] >= xx_arr[xx_off + 0]); // True
                                                                        // if
                                                                        // ascending
                                                                        // order
                                                                        // of
        // table, false otherwise.
        jl = 0; // Initialize lower
        ju = n - 1; // and upper limits.
        while (ju - jl > 1) { // If we are not yet done,
            jm = (ju + jl) >> 1; // compute a midpoint,
            if (x >= xx_arr[xx_off + jm] == ascnd)
                jl = jm; // and replace either the lower limit
            else
                ju = jm; // or the upper limit, as appropriate.
        } // Repeat until the test condition is satisfied.
        cor = abs(jl - jsav) > dj ? 0 : 1; // Decide whether to use hunt or
                                           // locate next time.
        jsav = jl;
        return MAX(0, MIN(n - mm, jl - ((mm - 2) >> 1)));
    }

    public int hunt(final double x) throws NRException {
        // Given a value x, return a value j such that x is (insofar as
        // possible) centered in the subrange xx[j..j+mm-1], where xx is the
        // stored pointer. The values in xx must be monotonic, either
        // increasing or decreasing. The returned value is not less than 0,
        // nor greater than n-1.

        int jl = jsav, jm, ju, inc = 1;
        if (n < 2 || mm < 2 || mm > n)
            throw new NRException("hunt size error");
        boolean ascnd = (xx_arr[xx_off + n - 1] >= xx_arr[xx_off + 0]); // True
                                                                        // if
                                                                        // ascending
                                                                        // order
                                                                        // of
        // table, false otherwise.
        if (jl < 0 || jl > n - 1) { // Input guess not useful. Go immediately to
                                    // bisec
            jl = 0; // tion.
            ju = n - 1;
        } else {
            if (x >= xx_arr[xx_off + jl] == ascnd) { // Hunt up:
                for (;;) {
                    ju = jl + inc;
                    if (ju >= n - 1) {
                        ju = n - 1;
                        break;
                    } // Off end of table.
                    else if (x < xx_arr[xx_off + ju] == ascnd)
                        break; // Found bracket.
                    else { // Not done, so double the increment and try again.
                        jl = ju;
                        inc += inc;
                    }
                }
            } else { // Hunt down:
                ju = jl;
                for (;;) {
                    jl = jl - inc;
                    if (jl <= 0) {
                        jl = 0;
                        break;
                    } // Off end of table.
                    else if (x >= xx_arr[xx_off + jl] == ascnd)
                        break; // Found bracket.
                    else { // Not done, so double the increment and try again.
                        ju = jl;
                        inc += inc;
                    }
                }
            }
        }
        while (ju - jl > 1) { // Hunt is done, so begin the final bisection
                              // phase:
            jm = (ju + jl) >> 1;
            if (x >= xx_arr[xx_off + jm] == ascnd)
                jl = jm;
            else
                ju = jm;
        }
        cor = abs(jl - jsav) > dj ? 0 : 1; // Decide whether to use hunt or
                                           // locate next
        jsav = jl; // time.
        return MAX(0, MIN(n - mm, jl - ((mm - 2) >> 1)));
    }

}
