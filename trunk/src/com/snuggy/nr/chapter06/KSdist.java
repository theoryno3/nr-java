package com.snuggy.nr.chapter06;

import static com.snuggy.nr.chapter06.Static.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class KSdist {

    // Kolmogorov-Smirnov cumulative distribution functions and their inverses.

    public double pks(final double z) throws NRException {
        // Return cumulative distribution function.
        if (z < 0.)
            throw new NRException("bad z in KSdist");
        if (z == 0.)
            return 0.;
        if (z < 1.18) {
            double y = exp(-1.23370055013616983 / SQR(z));
            return 2.25675833419102515 * sqrt(-log(y)) * (y + pow(y, 9) + pow(y, 25) + pow(y, 49));
        } else {
            double x = exp(-2. * SQR(z));
            return 1. - 2. * (x - pow(x, 4) + pow(x, 9));
        }
    }

    public double qks(final double z) throws NRException {
        // Return complementary cumulative distribution function.
        if (z < 0.)
            throw new NRException("bad z in KSdist");
        if (z == 0.)
            return 1.;
        if (z < 1.18)
            return 1. - pks(z);
        double x = exp(-2. * SQR(z));
        return 2. * (x - pow(x, 4) + pow(x, 9));
    }

    public double invqks(double q) throws NRException {
        // Return inverse of the complementary cumulative distribution function.
        double y, logy, x, xp, f, ff, u, t;
        @SuppressWarnings("unused")
        double yp;
        if (q <= 0. || q > 1.)
            throw new NRException("bad q in KSdist");
        if (q == 1.)
            return 0.;
        if (q > 0.3) {
            f = -0.392699081698724155 * SQR(1. - q);
            y = invxlogx(f); // Initial guess.
            do {
                yp = y;
                logy = log(y);
                ff = f / SQR(1. + pow(y, 4) + pow(y, 12));
                u = (y * logy - ff) / (1. + logy); // Newton’s method
                                                   // correction.
                y = y - (t = u / MAX(0.5, 1. - 0.5 * u / (y * (1. + logy)))); // Halley.
            } while (abs(t / y) > 1.e-15);
            return 1.57079632679489662 / sqrt(-log(y));
        } else {
            x = 0.03;
            do { // Iteration (6.14.59).
                xp = x;
                x = 0.5 * q + pow(x, 4) - pow(x, 9);
                if (x > 0.06)
                    x += pow(x, 16) - pow(x, 25);
            } while (abs((xp - x) / x) > 1.e-15);
            return sqrt(-0.5 * log(x));
        }
    }

    public double invpks(final double p) throws NRException {
        return invqks(1. - p);
    }
    // Return inverse of the cumulative distribution function.

}
