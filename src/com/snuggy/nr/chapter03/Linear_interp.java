
package com.snuggy.nr.chapter03;

public class Linear_interp extends Base_interp {

    // Piecewise linear interpolation object. Construct with x and y vectors,
    // then call interp for interpolated values.

    public Linear_interp(final double[] xv, final double[] yv) {
        super(xv, yv, 0, 2);
    }

    @Override
    public double rawinterp(final int j, final double x) {
        if (xx_arr[xx_off + j] == xx_arr[xx_off + j + 1])
            return yy_arr[yy_off + j]; // Table is defective, but we can
                                       // recover.
        else
            return yy_arr[yy_off + j] + ((x - xx_arr[xx_off + j]) / (xx_arr[xx_off + j + 1] - xx_arr[xx_off + j]))
                    * (yy_arr[yy_off + j + 1] - yy_arr[yy_off + j]);
    }

}
