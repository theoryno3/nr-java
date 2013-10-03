
package com.snuggy.nr.chapter17;

import com.snuggy.nr.util.*;

public interface Dtype {
    void eval(final double x, final double[] y, final double[] z) throws NRException;
    void jacobian(final double x, double[] y, double[] dfdx, double[][] dfdy) throws NRException;
}
