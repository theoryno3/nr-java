
package com.snuggy.nr.chapter03;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

public class RBF_inversemultiquadric extends RBF_fn {
    // Same as above, but for inverse multiquadric.
    private double r02;

    public RBF_inversemultiquadric() {
        this(1.0);
    }

    public RBF_inversemultiquadric(final double scale) {
        r02 = (SQR(scale));
    }

    @Override
    public double rbf(final double r) {
        return 1. / sqrt(SQR(r) + r02);
    }
}

