
package com.snuggy.nr.chapter06;

import static com.snuggy.nr.chapter06.Static.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Poissondist extends Gamma {

    // Poisson distribution, derived from the gamma function Gamma.
    private double lam;

    public Poissondist(final double llam) throws NRException {
        // Constructor. Initialize with .
        lam = (llam);
        if (lam <= 0.)
            throw new NRException("bad lam in Poissondist");
    }

    public double p(final int n) throws NRException {
        // Return probability density function.
        if (n < 0)
            throw new NRException("bad n in Poissondist");
        return exp(-lam + n * log(lam) - gammln(n + 1.));
    }

    public double cdf(final int n) throws NRException {
        // Return cumulative distribution function.
        if (n < 0)
            throw new NRException("bad n in Poissondist");
        if (n == 0)
            return 0.;
        return gammq((double) n, lam);
    }

    public int invcdf(final double p) throws NRException {
        // Given argument P, return integer n such that P.< n/  P  P.< nC1/.
        int n, nl, nu, inc = 1;
        if (p <= 0. || p >= 1.)
            throw new NRException("bad p in Poissondist");
        if (p < exp(-lam))
            return 0;
        n = (int) MAX(sqrt(lam), 5.); // Starting guess near peak of density.
        if (p < cdf(n)) { // Expand interval until we bracket.
            do {
                n = MAX(n - inc, 0);
                inc *= 2;
            } while (p < cdf(n));
            nl = n;
            nu = n + inc / 2;
        } else {
            do {
                n += inc;
                inc *= 2;
            } while (p > cdf(n));
            nu = n;
            nl = n - inc / 2;
        }
        while (nu - nl > 1) { // Now contract the interval by bisection.
            n = (nl + nu) / 2;
            if (p < cdf(n))
                nu = n;
            else
                nl = n;
        }
        return nl;
    }

}
