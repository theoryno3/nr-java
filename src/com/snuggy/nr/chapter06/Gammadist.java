
package com.snuggy.nr.chapter06;

import static com.snuggy.nr.chapter06.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Gammadist extends Gamma {

    // Gamma distribution, derived from the gamma function Gamma.
    private double alph, bet, fac;

    public Gammadist(final double aalph) throws NRException {
        this(aalph, 1.);
    }

    public Gammadist(final double aalph, final double bbet) throws NRException {
        // Constructor. Initialize with ? and ?.
        alph = (aalph);
        bet = (bbet);
        if (alph <= 0. || bet <= 0.)
            throw new NRException("bad alph,bet in Gammadist");
        fac = alph * log(bet) - gammln(alph);
    }

    public double p(final double x) throws NRException {
        // Return probability density function.
        if (x <= 0.)
            throw new NRException("bad x in Gammadist");
        return exp(-bet * x + (alph - 1.) * log(x) + fac);
    }

    public double cdf(final double x) throws NRException {
        // Return cumulative distribution function.
        if (x < 0.)
            throw new NRException("bad x in Gammadist");
        return gammp(alph, bet * x);
    }

    public double invcdf(final double p) throws NRException {
        // Return inverse cumulative distribution function.
        if (p < 0. || p >= 1.)
            throw new NRException("bad p in Gammadist");
        return invgammp(p, alph) / bet;
    }
}
