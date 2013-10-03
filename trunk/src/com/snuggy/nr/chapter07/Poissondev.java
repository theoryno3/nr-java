
package com.snuggy.nr.chapter07;

import static com.snuggy.nr.chapter06.Static.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Poissondev extends Ran {

    // Structure for Poisson deviates.
    private double lambda, sqlam, loglam, lamexp, lambold;
    private double[] logfact;
    @SuppressWarnings("unused")
    private int swch;

    public Poissondev(final double llambda, final long i) {
        super(i);
        lambda = (llambda);
        logfact = doub_arr(1024, -1.);
        lambold = (-1.);
    }

    // Constructor arguments are  and a random sequence seed.

    public int dev() throws NRException {
        // Return a Poisson deviate using the most recently set value of .
        double u, u2, v, v2 = 0.0, p, t, lfac;
        int k;
        if (lambda < 5.) { // Will use product of uniforms method.
            if (lambda != lambold)
                lamexp = exp(-lambda);
            k = -1;
            t = 1.;
            do {
                ++k;
                t *= doub();
            } while (t > lamexp);
        } else { // Will use ratio-of-uniforms method.
            if (lambda != lambold) {
                sqlam = sqrt(lambda);
                loglam = log(lambda);
            }
            for (;;) {
                u = 0.64 * doub();
                v = -0.68 + 1.28 * doub();
                if (lambda > 13.5) { // Outer squeeze for fast rejection.
                    v2 = SQR(v);
                    if (v >= 0.) {
                        if (v2 > 6.5 * u * (0.64 - u) * (u + 0.2))
                            continue;
                    } else {
                        if (v2 > 9.6 * u * (0.66 - u) * (u + 0.07))
                            continue;
                    }
                }
                k = Int(floor(sqlam * (v / u) + lambda + 0.5));
                if (k < 0)
                    continue;
                u2 = SQR(u);
                if (lambda > 13.5) { // Inner squeeze for fast acceptance.
                    if (v >= 0.) {
                        if (v2 < 15.2 * u2 * (0.61 - u) * (0.8 - u))
                            break;
                    } else {
                        if (v2 < 6.76 * u2 * (0.62 - u) * (1.4 - u))
                            break;
                    }
                }
                if (k < 1024) {
                    if (logfact[k] < 0.)
                        logfact[k] = gammln(k + 1.);
                    lfac = logfact[k];
                } else
                    lfac = gammln(k + 1.);
                p = sqlam * exp(-lambda + k * loglam - lfac); // Only when we
                                                              // must.
                if (u2 < p)
                    break;
            }
        }
        lambold = lambda;
        return k;
    }

    public int dev(final double llambda) throws NRException {
        // Reset  and then return a Poisson deviate.
        lambda = llambda;
        return dev();
    }

}
