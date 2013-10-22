
package com.snuggy.nr.chapter09;

import static com.snuggy.nr.util.Static.*;
import static com.snuggy.nr.refs.Refs.*;

import com.snuggy.nr.refs.*;
import com.snuggy.nr.util.*;

public class NRfmin<T extends Func_DoubVec_To_DoubVec> implements Func_DoubVec_To_Doub {

    // Returns f D 1
    // 2F  F. Also stores value of F in fvec.
    private final T func;
    private $double1d fvec;
    private int n;

    public NRfmin(final T funcc) throws NRException {
        // Initialize with user-supplied function or functor that returns the
        // vector of functions to be zeroed.
        func = (funcc);
        fvec = $(doub_vec(0));
    }

    public double eval(final double[] x) throws NRException {
        // Returns f at x, and stores F.x/ in fvec.
        n = x.length;
        double sum = 0;
        fvec.$(func.eval(x));
        for (int i = 0; i < n; i++)
            sum += SQR(fvec.$()[i]);
        return 0.5 * sum;
    }
    
    public $double1d fvec() {
        return fvec;
    }

}
