
package com.snuggy.nr.chapter17;

import com.snuggy.nr.util.*;

public class StepperStoerm_rhs implements Dtype {
    @Override
    public void eval(final double x, final double[] y, final double[] dydx) {
        dydx[0] = x * y[0];
        dydx[1] = x * y[1];
    }
    @Override
    public void jacobian(final double x, final double[] y, final double[] dfdx, final double[][] dfdy) throws NRException {
        throw new NRException();
    }
}