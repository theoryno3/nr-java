
package com.snuggy.nr.chapter17;

import com.snuggy.nr.util.*;

// TODO: can remove this?
public interface IStepperBS extends IStepperBase {
    void step(double htry, Dtype derivs) throws NRException;

    boolean dy(double[] y, double htot, int k, double[] yend, int ipt_ref[], Dtype derivs)
            throws NRException;

    void polyextr(int k, double[][] table, double[] last);

    void prepare_dense(double h, double[] dydxnew, double[] ysav, double[] scale, int k,
            double error_ref[]);

    double dense_out(int i, double x, double h) throws NRException;

    void dense_interp(int n, double[] y, int imit);
    
}