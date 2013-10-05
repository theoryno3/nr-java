
package com.snuggy.nr.chapter17;

import com.snuggy.nr.refs.*;
import com.snuggy.nr.util.*;

// TODO: can remove this?
public interface IStepperBS extends IStepperBase {
    void step(double htry, Dtype derivs) throws NRException;

    boolean dy(final double[] y, double htot, int k, final double[] yend, int ipt_ref[], Dtype derivs)
            throws NRException;

    void polyextr(int k, final double[][] table, final double[] last);

    void prepare_dense(double h, final double[] dydxnew, final double[] ysav, final double[] scale, int k,
            $double error);

    double dense_out(int i, double x, double h) throws NRException;

    void dense_interp(int n, final double[] y, int imit);
    
}