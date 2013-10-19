
package com.snuggy.nr.chapter03;

import static com.snuggy.nr.refs.Refs.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.refs.*;
import com.snuggy.nr.util.*;

public class Poly_interp extends Base_interp {

    // Polynomial interpolation object. Construct with x and y vectors,
    // and the number M of points to be used locally (polynomial order plus
    // one), then call interp for interpolated values.

    protected double dy;
    
    public double dy() {
        return dy;
    }

    public Poly_interp(final double[] xv, final double[] yv, final int m) {
        super(xv, $_(yv, 0), m);
        dy = (0.0);
    }

    // public double rawinterp(int jl, double x)

    public double rawinterp(final int jl, final double x) throws NRException {
        // Given a value x, and using pointers to data xx and yy, this routine
        // returns an interpolated value y, and stores an error estimate dy.
        // The returned value is obtained by mm-point polynomial interpolation
        // on the subrange xx[jl..jl+mm-1].

        int i, m, ns = 0;
        double y, den, dif, dift, ho, hp, w;
        final $double xa = $(xx, jl); 
        final $double ya = $(yy, jl);
        final double[] c = doub_vec(mm), d = doub_vec(mm);
        dif = abs(x - xa.$_(0));
        for (i = 0; i < mm; i++) { // Here we find the index ns of the closest
                                   // table entry,
            if ((dift = abs(x - xa.$_(i))) < dif) {
                ns = i;
                dif = dift;
            }
            c[i] = ya.$_(i); // and initialize the tableau of c’s and
                                       // d’s.
            d[i] = ya.$_(i);
        }
        y = ya.$_((ns--)); // This is the initial approximation to y.
        for (m = 1; m < mm; m++) { // For each column of the tableau,
            for (i = 0; i < mm - m; i++) { // we loop over the current c’s and
                                           // d’s and update
                ho = xa.$_(i) - x; // them.
                hp = xa.$_(i + m) - x;
                w = c[i + 1] - d[i];
                if ((den = ho - hp) == 0.0)
                    throw new NRException("Poly_interp error");
                // This error can occur only if two input xa’s are (to within
                // roundoff) identical.
                den = w / den;
                d[i] = hp * den; // Here the c’s and d’s are updated.
                c[i] = ho * den;
            }
            y += (dy = (2 * (ns + 1) < (mm - m) ? c[ns + 1] : d[ns--]));
            // After each column in the tableau is completed, we decide which
            // correction, c or d, we want to add to our accumulating value of
            // y, i.e., which path to take through the tableau — forking up or
            // down. We do this in such a way as to take the most “straight
            // line”
            // route through the tableau to its apex, updating ns accordingly to
            // keep track of where we are. This route keeps the partial
            // approximations centered (insofar as possible) on the target x.
            // The last dy added is thus the error indication.
        }
        return y;
    }

}
