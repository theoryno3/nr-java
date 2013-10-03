
package com.snuggy.nr.chapter06;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Lognormaldist extends Erf {

    // Lognormal distribution, derived from the error function Erf.
    private double mu, sig;

    public Lognormaldist() throws NRException {
        this(0., 1.);
    }

    public Lognormaldist(final double mmu) throws NRException {
        this(mmu, 1.);
    }

    public Lognormaldist(final double mmu, final double ssig) throws NRException {
        mu = (mmu);
        sig = (ssig);
        if (sig <= 0.)
            throw new NRException("bad sig in Lognormaldist");
    }

    public double p(final double x) throws NRException {
        // Return probability density function.
        if (x < 0.)
            throw new NRException("bad x in Lognormaldist");
        if (x == 0.)
            return 0.;
        return (0.398942280401432678 / (sig * x)) * exp(-0.5 * SQR((log(x) - mu) / sig));
    }

    public double cdf(final double x) throws NRException {
        // Return cumulative distribution function.
        if (x < 0.)
            throw new NRException("bad x in Lognormaldist");
        if (x == 0.)
            return 0.;
        return 0.5 * erfc(-0.707106781186547524 * (log(x) - mu) / sig);
    }

    public double invcdf(final double p) throws NRException {
        // Return inverse cumulative distribution function.
        if (p <= 0. || p >= 1.)
            throw new NRException("bad p in Lognormaldist");
        return exp(-1.41421356237309505 * sig * inverfc(2. * p) + mu);
    }

}
