
package com.snuggy.nr.chapter17;

import com.snuggy.nr.util.*;

public class StepperBS_rhs implements Dtype {
    @Override
    public void eval(final double x, final double[] y, final double[] dydx) {
        dydx[0] = -y[1];
        dydx[1] = y[0] - (1.0 / x) * y[1];
        dydx[2] = y[1] - (2.0 / x) * y[2];
        dydx[3] = y[2] - (3.0 / x) * y[3];
    }
    @Override
    public void jacobian(final double x, final double[] y, final double[] dfdx, final double[][] dfdy) throws NRException {
        throw new NRException();
    }
}