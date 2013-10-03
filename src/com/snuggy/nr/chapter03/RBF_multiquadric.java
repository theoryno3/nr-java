
package com.snuggy.nr.chapter03;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

public class RBF_multiquadric extends RBF_fn {

    @Override
    public double rbf(final double r) {
        return sqrt(SQR(r) + r02);
    }

    // Instantiate this and send to RBF_interp to get multiquadric
    // interpolation.
    private double r02;

    public RBF_multiquadric() {
        this(1.0);
    }

    public RBF_multiquadric(final double scale) {
        r02 = (SQR(scale));
    }

    // Constructor argument is the scale factor. See text.

}
