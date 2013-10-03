
package com.snuggy.nr.chapter07;

import static com.snuggy.nr.util.Static.*;

public class Cauchydev extends Ran {

    // Structure for Cauchy deviates.
    private double mu, sig;

    public Cauchydev(final double mmu, final double ssig, final long i) {
        super(i);
        mu = (mmu);
        sig = (ssig);
    }

    // Constructor arguments are , , and a random sequence seed.

    public double dev() {
        // Return a Cauchy deviate.
        double v1, v2;
        do { // Find a random point in the unit semicircle.
            v1 = 2.0 * doub() - 1.0;
            v2 = doub();
        } while (SQR(v1) + SQR(v2) >= 1. || v2 == 0.);
        return mu + sig * v1 / v2; // Ratio of its coordinates is the tangent of
                                   // a
    } // random angle.

}
