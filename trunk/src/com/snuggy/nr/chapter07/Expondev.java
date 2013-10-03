
package com.snuggy.nr.chapter07;

import static java.lang.Math.*;

public class Expondev extends Ran {

    // Structure for exponential deviates.
    private double beta;

    public Expondev(final double bbeta, long i) {
        super(i);
        beta = (bbeta);
    }

    // Constructor arguments are ? and a random sequence seed.

    public double dev() {
        // Return an exponential deviate.
        double u;
        do
            u = doub();
        while (u == 0.);
        return -log(u) / beta;
    }

}
