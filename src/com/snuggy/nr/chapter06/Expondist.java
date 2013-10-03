
package com.snuggy.nr.chapter06;

import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Expondist {

    // Exponential distribution.
    private double bet;

    public Expondist(final double bbet) throws NRException {
        // Constructor. Initialize with ?.
        bet = (bbet);
        if (bet <= 0.)
            throw new NRException("bad bet in Expondist");
    }

    public double p(final double x) throws NRException {
        // Return probability density function.
        if (x < 0.)
            throw new NRException("bad x in Expondist");
        return bet * exp(-bet * x);
    }

    public double cdf(final double x) throws NRException {
        // Return cumulative distribution function.
        if (x < 0.)
            throw new NRException("bad x in Expondist");
        return 1. - exp(-bet * x);
    }

    public double invcdf(final double p) throws NRException {
        // Return inverse cumulative distribution function.
        if (p < 0. || p >= 1.)
            throw new NRException("bad p in Expondist");
        return -log(1. - p) / bet;
    }

}
