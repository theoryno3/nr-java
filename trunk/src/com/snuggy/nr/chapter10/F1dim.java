
package com.snuggy.nr.chapter10;

import static com.snuggy.nr.util.Static.*;

import com.snuggy.nr.util.*;

public class F1dim<T extends Func_DoubVec_To_Doub> implements Func_Doub_To_Doub {

    // Must accompany linmin in Linemethod.
    private final double[] p;
    private final double[] xi;
    private final int n;
    private final T func;
    private final double[] xt;

    public F1dim(final double[] pp, final double[] xii, final T funcc) {
        // Constructor takes as inputs an n-dimensional point p[0..n-1] and
        // an n-dimensional direction xi[0..n-1] from linmin, as well as the
        // function or functor that takes a vector argument.
        p = (pp);
        xi = (xii);
        n = (pp.length);
        func = (funcc);
        xt = doub_vec(n);
    }

    public double eval(final double x) throws NRException {
        // Functor returning value of the given function along a one-dimensional
        // line.
        for (int j = 0; j < n; j++)
            xt[j] = p[j] + x * xi[j];
        return func.eval(xt);
    }

}
