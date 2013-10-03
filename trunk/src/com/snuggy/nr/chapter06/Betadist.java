
package com.snuggy.nr.chapter06;

import static com.snuggy.nr.chapter06.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Betadist extends Beta {

    // Beta distribution, derived from the beta function Beta.
    private double alph, bet, fac;

    public Betadist(final double aalph, final double bbet) throws NRException {
        // Constructor. Initialize with ? and ?.
        alph = (aalph);
        bet = (bbet);
        if (alph <= 0. || bet <= 0.)
            throw new NRException("bad alph,bet in Betadist");
        fac = gammln(alph + bet) - gammln(alph) - gammln(bet);
    }

    public double p(final double x) throws NRException {
        // Return probability density function.
        if (x <= 0. || x >= 1.)
            throw new NRException("bad x in Betadist");
        return exp((alph - 1.) * log(x) + (bet - 1.) * log(1. - x) + fac);
    }

    public double cdf(final double x) throws NRException {
        // Return cumulative distribution function.
        if (x < 0. || x > 1.)
            throw new NRException("bad x in Betadist");
        return betai(alph, bet, x);
    }

    public double invcdf(final double p) throws NRException {
        // Return inverse cumulative distribution function.
        if (p < 0. || p > 1.)
            throw new NRException("bad p in Betadist");
        return invbetai(p, alph, bet);
    }

}
