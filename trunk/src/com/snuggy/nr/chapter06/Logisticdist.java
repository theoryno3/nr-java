
package com.snuggy.nr.chapter06;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Logisticdist {

    // Logistic distribution.
    private double mu, sig;

    public Logisticdist() throws NRException {
        this(0., 1.);
    }

    public Logisticdist(final double mmu) throws NRException {
        this(mmu, 1.);
    }

    public Logisticdist(final double mmu, final double ssig) throws NRException {
        // Constructor. Initialize with  and . The default with no arguments
        // is Logistic.0; 1/.
        mu = (mmu);
        sig = (ssig);
        if (sig <= 0.)
            throw new NRException("bad sig in Logisticdist");
    }

    public double p(final double x) {
        // Return probability density function.
        double e = exp(-abs(1.81379936423421785 * (x - mu) / sig));
        return 1.81379936423421785 * e / (sig * SQR(1. + e));
    }

    public double cdf(final double x) {
        // Return cumulative distribution function.
        double e = exp(-abs(1.81379936423421785 * (x - mu) / sig));
        if (x >= mu)
            return 1. / (1. + e); // Because we used abs to control overelse
        return e / (1. + e); // flow, we now have two cases.
    }

    public double invcdf(final double p) throws NRException {
        // Return inverse cumulative distribution function.
        if (p <= 0. || p >= 1.)
            throw new NRException("bad p in Logisticdist");
        return mu + 0.551328895421792049 * sig * log(p / (1. - p));
    }

}
