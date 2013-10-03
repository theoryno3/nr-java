
package com.snuggy.nr.chapter06;

import static com.snuggy.nr.chapter06.Static.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Studenttdist extends Beta {

    // Student-t distribution, derived from the beta function Beta.
    private double nu, mu, sig, np, fac;

    public Studenttdist(final double nnu) throws NRException {
        this(nnu, 0., 1.);
    }

    public Studenttdist(final double nnu, final double mmu) throws NRException {
        this(nnu, mmu, 1.);
    }

    public Studenttdist(final double nnu, final double mmu, final double ssig) throws NRException {
        // Constructor. Initialize with ,  and . The default with one
        // argument
        // is Student.; 0; 1/.
        nu = (nnu);
        mu = (mmu);
        sig = (ssig);
        if (sig <= 0. || nu <= 0.)
            throw new NRException("bad sig,nu in Studentdist");
        np = 0.5 * (nu + 1.);
        fac = gammln(np) - gammln(0.5 * nu);
    }

    public double p(final double t) {
        // Return probability density function.
        return exp(-np * log(1. + SQR((t - mu) / sig) / nu) + fac) / (sqrt(3.14159265358979324 * nu) * sig);
    }

    public double cdf(final double t) throws NRException {
        // Return cumulative distribution function.
        double p = 0.5 * betai(0.5 * nu, 0.5, nu / (nu + SQR((t - mu) / sig)));
        if (t >= mu)
            return 1. - p;
        else
            return p;
    }

    public double invcdf(final double p) throws NRException {
        // Return inverse cumulative distribution function.
        if (p <= 0. || p >= 1.)
            throw new NRException("bad p in Studentdist");
        double x = invbetai(2. * MIN(p, 1. - p), 0.5 * nu, 0.5);
        x = sig * sqrt(nu * (1. - x) / x);
        return (p >= 0.5 ? mu + x : mu - x);
    }

    public double aa(final double t) throws NRException {
        // Return the two-tailed cdf A.tj/.
        if (t < 0.)
            throw new NRException("bad t in Studentdist");
        return 1. - betai(0.5 * nu, 0.5, nu / (nu + SQR(t)));
    }

    public double invaa(final double p) throws NRException {
        // Return the inverse, namely t such that p D A.tj/.
        if (p < 0. || p >= 1.)
            throw new NRException("bad p in Studentdist");
        double x = invbetai(1. - p, 0.5 * nu, 0.5);
        return sqrt(nu * (1. - x) / x);
    }

}
