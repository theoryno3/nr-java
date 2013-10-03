
package com.snuggy.nr.chapter06;

import static com.snuggy.nr.chapter06.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Fdist extends Beta {

    // F distribution, derived from the beta function Beta.
    private double nu1, nu2;
    private double fac;

    public Fdist(final double nnu1, final double nnu2) throws NRException {
        // Constructor. Initialize with 1 and 2.
        nu1 = (nnu1);
        nu2 = (nnu2);
        if (nu1 <= 0. || nu2 <= 0.)
            throw new NRException("bad nu1,nu2 in Fdist");
        fac = 0.5 * (nu1 * log(nu1) + nu2 * log(nu2)) + gammln(0.5 * (nu1 + nu2)) - gammln(0.5 * nu1)
                - gammln(0.5 * nu2);
    }

    public double p(final double f) throws NRException {
        // Return probability density function.
        if (f <= 0.)
            throw new NRException("bad f in Fdist");
        return exp((0.5 * nu1 - 1.) * log(f) - 0.5 * (nu1 + nu2) * log(nu2 + nu1 * f) + fac);
    }

    public double cdf(final double f) throws NRException {
        // Return cumulative distribution function.
        if (f < 0.)
            throw new NRException("bad f in Fdist");
        return betai(0.5 * nu1, 0.5 * nu2, nu1 * f / (nu2 + nu1 * f));
    }

    public double invcdf(final double p) throws NRException {
        // Return inverse cumulative distribution function.
        if (p <= 0. || p >= 1.)
            throw new NRException("bad p in Fdist");
        double x = invbetai(p, 0.5 * nu1, 0.5 * nu2);
        return nu2 * x / (nu1 * (1. - x));
    }

}
