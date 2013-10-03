
package com.snuggy.nr.chapter06;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Cauchydist {

    // Cauchy distribution.
    private double mu, sig;

    public Cauchydist() throws NRException {
        this(0., 1.);
    }

    public Cauchydist(final double mmu) throws NRException {
        this(mmu, 1.);
    }

    public Cauchydist(final double mmu, final double ssig) throws NRException {
        // Constructor. Initialize with  and . The default with no arguments
        // is Cauchy.0; 1/.
        mu = (mmu);
        sig = (ssig);
        if (sig <= 0.)
            throw new NRException("bad sig in Cauchydist");
    }

    public double p(final double x) {
        // Return probability density function.
        return 0.318309886183790671 / (sig * (1. + SQR((x - mu) / sig)));
    }

    public double cdf(final double x) {
        // Return cumulative distribution function.
        return 0.5 + 0.318309886183790671 * atan2(x - mu, sig);
    }

    public double invcdf(final double p) throws NRException {
        // Return inverse cumulative distribution function.
        if (p <= 0. || p >= 1.)
            throw new NRException("bad p in Cauchydist");
        return mu + sig * tan(3.14159265358979324 * (p - 0.5));
    }

}
