
package com.snuggy.nr.chapter07;

import static java.lang.Math.*;

public class Normaldev_BM extends Ran {

    // Structure for normal deviates.
    private double mu, sig;
    private double storedval;

    public Normaldev_BM(final double mmu, final double ssig, final long i) {
        super(i);
        mu = (mmu);
        sig = (ssig);
        storedval = (0.);
    }

    // Constructor arguments are , , and a random sequence seed.

    public double dev() {
        // Return a normal deviate.
        double v1, v2, rsq, fac;
        if (storedval == 0.) { // We don’t have an extra deviate handy, so
            do {
                v1 = 2.0 * doub() - 1.0; // pick two uniform numbers in the
                                         // square ex
                v2 = 2.0 * doub() - 1.0; // tending from -1 to +1 in each
                                         // direction,
                rsq = v1 * v1 + v2 * v2; // see if they are in the unit circle,
            } while (rsq >= 1.0 || rsq == 0.0); // or try again.
            fac = sqrt(-2.0 * log(rsq) / rsq); // Now make the Box-Muller
                                               // transformation to
            // get two normal deviates. Return one and
            // save the other for next time.
            storedval = v1 * fac;
            return mu + sig * v2 * fac;
        } else { // We have an extra deviate handy,
            fac = storedval;
            storedval = 0.;
            return mu + sig * fac; // so return it.
        }
    }

}
