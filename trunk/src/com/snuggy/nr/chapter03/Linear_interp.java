
package com.snuggy.nr.chapter03;

import static com.snuggy.nr.refs.Refs.*;

import com.snuggy.nr.util.*;

public class Linear_interp extends Base_interp {

    // Piecewise linear interpolation object. Construct with x and y vectors,
    // then call interp for interpolated values.

    public Linear_interp(final double[] xv, final double[] yv) {
        super(xv, $(yv, 0), 2);
    }

    @Override
    public double rawinterp(final int j, final double x) throws NRException {
        if (xx.$(j) == xx.$(j + 1))
            return yy.$(j); // Table is defective, but we can
                                       // recover.
        else
            return yy.$(j) + ((x - xx.$(j)) / (xx.$(j + 1) - xx.$(j)))
                    * (yy.$(j + 1) - yy.$(j));
    }

}
