
package com.snuggy.nr.chapter07;

import static java.lang.Math.*;

public class Logisticdev extends Ran {

    // Structure for logistic deviates.
    private double mu, sig;

    public Logisticdev(final double mmu, final double ssig, final long i) {
        super(i);
        mu = (mmu);
        sig = (ssig);
    }

    // Constructor arguments are , , and a random sequence seed.

    public double dev() {
        // Return a logistic deviate.
        double u;
        do
            u = doub();
        while (u * (1. - u) == 0.);
        return mu + 0.551328895421792050 * sig * log(u / (1. - u));
    }

}
