
package com.snuggy.nr.chapter15;

import static com.snuggy.nr.refs.Refs.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.refs.*;
import com.snuggy.nr.util.*;

public class Static {

    private static int fpoly_np = 10; // Global variable for the degree plus
                                      // one. fit examples.h

    public static final double[] fpoly(final double x) {
        // Fitting routine for a polynomial of degree fpoly_np-1.
        int j;
        final double[] p = doub_vec(fpoly_np);
        p[0] = 1.0;
        for (j = 1; j < fpoly_np; j++)
            p[j] = p[j - 1] * x;
        return p;
    }

    // The second example is slightly less trivial. It is used to fit Legendre
    // polynomials up to some order fleg_nl to a data set. (Note that, for most
    // uses, the data should satisfy 1 x 1.)

    private static int fleg_nl = 10; // Global variable for the degree plus one.
                                     // fit examples.h

    public static final double[] fleg(final double x) {
        // Fitting routine for an expansion with nl Legendre polynomials,
        // evaluated using the recurrence relation as in 5.4.
        int j;
        double twox, f2, f1, d;
        final double[] pl = doub_vec(fleg_nl);
        pl[0] = 1.;
        pl[1] = x;
        if (fleg_nl > 2) {
            twox = 2. * x;
            f2 = x;
            d = 1.;
            for (j = 2; j < fleg_nl; j++) {
                f1 = d++;
                f2 += twox;
                pl[j] = (f2 * pl[j - 1] - f1 * pl[j - 2]) / d;
            }
        }
        return pl;
    }

    public static final double[] quadratic2d(final double[] xx) {
        final double[] ans = doub_vec(6);
        double x = xx[0], y = xx[1];
        ans[0] = 1;
        ans[1] = x;
        ans[2] = y;
        ans[3] = x * x;
        ans[4] = x * y;
        ans[5] = y * y;
        return ans;
    }

    public static void fgauss(final double x, final double[] a, 
            final $double y, final double[] dyda) {
        // y.xI a/ is the sum of na/3 Gaussians (15.5.16). The amplitude,
        // center, and width of the Gaussians are stored in consecutive
        // locations of a: a[3k] D Bk, a[3kC1] D Ek, a[3kC2] D Gk, k D 0; :::;
        // na/3  1. The dimensions of the arrays are a[0..na-1], dyda[0..na-1].
        int i, na = a.length;
        double fac, ex, arg;
        y.$(0.);
        for (i = 0; i < na - 1; i += 3) {
            arg = (x - a[i + 1]) / a[i + 2];
            ex = exp(-SQR(arg));
            fac = a[i] * ex * 2. * arg;
            y.$(y.$() + a[i] * ex);
            dyda[i] = ex;
            dyda[i + 1] = fac / a[i + 2];
            dyda[i + 2] = fac * arg / a[i + 2];
        }
    }

    public static double mcmcstep(final int m, final $$<State> s, final Plog plog, final Proposal propose)
            throws NRException {
        // Take m MCMC steps, starting with (and updating) s.
        State sprop = new State(); // Storage for candidate.
        double alph, ran;
        $double qratio = $(0.0);
        int accept = 0;
        plog.func(s.$());
        for (int i = 0; i < m; i++) { // Loop over steps.
            propose.func(s.$(), sprop, qratio);
            alph = min(1., qratio.$() * exp(plog.func(sprop) - s.$().plog())); // Equation
                                                                                     // (15.8.5).
            ran = propose.gau().doub();
            if (ran < alph) { // Accept the candidate.
                $$(s, sprop);
                plog.func(s.$());
                accept++;
            }
        }
        return accept / Doub(m);
    }
}
