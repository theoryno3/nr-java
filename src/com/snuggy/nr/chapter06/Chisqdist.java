
package com.snuggy.nr.chapter06;

import static com.snuggy.nr.chapter06.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Chisqdist extends Gamma {

    // 2 distribution, derived from the gamma function Gamma.
    private double nu, fac;

    public Chisqdist(final double nnu) throws NRException {
        // Constructor. Initialize with .
        nu = (nnu);
        if (nu <= 0.)
            throw new NRException("bad nu in Chisqdist");
        fac = 0.693147180559945309 * (0.5 * nu) + gammln(0.5 * nu);
    }

    public double p(final double x2) throws NRException {
        // Return probability density function.
        if (x2 <= 0.)
            throw new NRException("bad x2 in Chisqdist");
        return exp(-0.5 * (x2 - (nu - 2.) * log(x2)) - fac);
    }

    public double cdf(final double x2) throws NRException {
        // Return cumulative distribution function.
        if (x2 < 0.)
            throw new NRException("bad x2 in Chisqdist");
        return gammp(0.5 * nu, 0.5 * x2);
    }

    public double invcdf(final double p) throws NRException {
        // Return inverse cumulative distribution function.
        if (p < 0. || p >= 1.)
            throw new NRException("bad p in Chisqdist");
        return 2. * invgammp(p, 0.5 * nu);
    }

}
