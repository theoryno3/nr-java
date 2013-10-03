
package com.snuggy.nr.chapter03;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

public class RBF_thinplate extends RBF_fn {
    // Same as above, but for thin-plate spline.
    private double r0;

    public RBF_thinplate() {
        this(1.0);
    }

    public RBF_thinplate(final double scale) {
        r0 = (scale);
    }

    @Override
    public double rbf(final double r) {
        return r <= 0. ? 0. : SQR(r) * log(r / r0);
    }
}

