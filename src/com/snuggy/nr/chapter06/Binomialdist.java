
package com.snuggy.nr.chapter06;

import static com.snuggy.nr.chapter06.Static.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Binomialdist extends Beta {

    // Binomial distribution, derived from the beta function Beta.
    private int n;
    private double pe, fac;

    public Binomialdist(final int nn, final double ppe) throws NRException {
        // Constructor. Initialize with n (sample size) and p (event
        // probability).
        n = (nn);
        pe = (ppe);
        if (n <= 0 || pe <= 0. || pe >= 1.)
            throw new NRException("bad args in Binomialdist");
        fac = gammln(n + 1.);
    }

    public double p(final int k) throws NRException {
        // Return probability density function.
        if (k < 0)
            throw new NRException("bad k in Binomialdist");
        if (k > n)
            return 0.;
        return exp(k * log(pe) + (n - k) * log(1. - pe) + fac - 
                   gammln(k + 1.) - gammln(n - k + 1.));
    }

    public double cdf(final int k) throws NRException {
        // Return cumulative distribution function.
        if (k < 0)
            throw new NRException("bad k in Binomialdist");
        if (k == 0)
            return 0.;
        if (k > n)
            return 1.;
        return 1. - betai((double) k, n - k + 1., pe);
    }

    public int invcdf(final double p) throws NRException {
        // Given argument P, return integer n such that P.< n/  P  P.< nC1/.
        int k, kl, ku, inc = 1;
        if (p <= 0. || p >= 1.)
            throw new NRException("bad p in Binomialdist");
        k = MAX(0, MIN(n, (int) (n * pe))); // Starting guess near peak of
                                            // density.
        if (p < cdf(k)) { // Expand interval until we bracket.
            do {
                k = MAX(k - inc, 0);
                inc *= 2;
            } while (p < cdf(k));
            kl = k;
            ku = k + inc / 2;
        } else {
            do {
                k = MIN(k + inc, n + 1);
                inc *= 2;
            } while (p > cdf(k));
            ku = k;
            kl = k - inc / 2;
        }
        while (ku - kl > 1) { // Now contract the interval by bisection.
            k = (kl + ku) / 2;
            if (p < cdf(k))
                ku = k;
            else
                kl = k;
        }
        return kl;
    }

}
