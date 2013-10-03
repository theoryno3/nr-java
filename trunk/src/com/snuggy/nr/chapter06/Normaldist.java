
package com.snuggy.nr.chapter06;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Normaldist extends Erf {

    // Normal distribution, derived from the error function Erf.
    private double mu, sig;

    public Normaldist() throws NRException {
        this(0., 1.);
    }

    public Normaldist(final double mmu) throws NRException {
        this(mmu, 1.);
    }

    public Normaldist(final double mmu, final double ssig) throws NRException {
        mu = (mmu);
        sig = (ssig);
        // Constructor. Initialize with  and . The default with
        // no arguments is N.0; 1/.
        if (sig <= 0.)
            throw new NRException("bad sig in Normaldist");
    }

    public double p(final double x) {
        // Return probability density function.
        return (0.398942280401432678 / sig) * exp(-0.5 * SQR((x - mu) / sig));
    }

    public double cdf(final double x) throws NRException {
        // Return cumulative distribution function.
        return 0.5 * erfc(-0.707106781186547524 * (x - mu) / sig);
    }

    public double invcdf(final double p) throws NRException {
        // Return inverse cumulative distribution function.
        if (p <= 0. || p >= 1.)
            throw new NRException("bad p in Normaldist");
        return -1.41421356237309505 * sig * inverfc(2. * p) + mu;
    }

}
