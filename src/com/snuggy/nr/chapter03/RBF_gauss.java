
package com.snuggy.nr.chapter03;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

public class RBF_gauss extends RBF_fn {
    // Same as above, but for Gaussian.
    private double r0;

    public RBF_gauss() {
        this(1.0);
    }

    public RBF_gauss(final double scale) {
        r0 = (scale);
    }

    @Override
    public double rbf(final double r) {
        return exp(-0.5 * SQR(r / r0));
    }
}

