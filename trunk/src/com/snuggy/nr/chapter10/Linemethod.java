
package com.snuggy.nr.chapter10;

import static com.snuggy.nr.util.Static.*;

import com.snuggy.nr.util.*;

import static com.snuggy.nr.refs.Refs.*;

import com.snuggy.nr.refs.*;

public class Linemethod<T extends Func_DoubVec_To_Doub> {

    // Base class for line-minimization algorithms. Provides the
    // line-minimization
    // routine linmin.
    private final T func;
    private $$double1d p;
    private $$double1d xi;
    private int n;

    public Linemethod(final T funcc) throws NRException {
        // Constructor argument is the user-supplied function or functor to be
        // minimized.
        p = $$(doub_vec(0));
        xi = $$(doub_vec(0));
        func = (funcc);
    }
    
    public double[] p() {
        return p.$();
    }
    
    public double[] xi() {
        return xi.$();
    }
    
    public void p(double[] p) throws NRException {
        $$(this.p, p);
    }
    
    public void xi(double[] xi) throws NRException {
        $$(this.xi, xi);
    }

    public double linmin() throws NRException {
        // Line-minimization routine. Given an n-dimensional point p[0..n-1]
        // and an n-dimensional direction xi[0..n-1], moves and resets p to
        // where
        // the function or functor func(p) takes on a minimum along the
        // direction
        // xi from p, and replaces xi by the actual vector displacement that p
        // was moved. Also returns the value of func at the returned location p.
        // This is actually all accomplished by calling the routines bracket
        // and minimize of Brent.
        double ax, xx, xmin;
        n = p.$().length;
        F1dim<T> f1dim = new F1dim<T>(p.$(), xi.$(), func);
        ax = 0.0; // Initial guess for brackets.
        xx = 1.0;
        Brent brent = new Brent();
        brent.bracket(ax, xx, f1dim);
        xmin = brent.minimize(f1dim);
        for (int j = 0; j < n; j++) { // Construct the vector results to return.
            xi.$()[j] *= xmin;
            p.$()[j] += xi.$()[j];
        }
        return brent.fmin();
    }

}
