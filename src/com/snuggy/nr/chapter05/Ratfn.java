
package com.snuggy.nr.chapter05;

import static com.snuggy.nr.refs.Refs.*;
import static com.snuggy.nr.util.Static.*;

import com.snuggy.nr.refs.*;
import com.snuggy.nr.util.*;

public class Ratfn implements ByValue<Ratfn> {

    // Function object for a rational function.

    private $$double1d cofs;
    private int nn, dd; // Number of numerator, denominator coefficients.
    
    private Ratfn() {
    }

    @Override
    public Ratfn copyOut() throws NRException {
        Ratfn r = new Ratfn();
        r.cofs = $$(cofs.$());
	    r.nn = nn;
	    r.dd = dd; 
	    return r;
    }

    @Override
    public void copyIn(Ratfn t) throws NRException {
        $$(cofs, t.cofs);
	    nn = t.nn;
	    dd = t.dd; 
    }

    public Ratfn(final double[] num, final double[] den) throws NRException {
        // Constructor from numerator, denominator polyomials (as coefficient
        // vectors).
        cofs = $$(doub_vec(num.length + den.length - 1));
        nn = (num.length);
        dd = (den.length);
        int j;
        for (j = 0; j < nn; j++)
            cofs.$()[j] = num[j] / den[0];
        for (j = 1; j < dd; j++)
            cofs.$()[j + nn - 1] = den[j] / den[0];
    }

    public Ratfn(final double[] coffs, final int n, final int d) throws NRException {
        cofs = $$(coffs);
        nn = (n);
        dd = (d);
    } // Constructor from coefficients already normalized and in a single array.

    public double func(final double x) {
        // Evaluate the rational function at x and return result.
        int j;
        double sumn = 0., sumd = 0.;
        for (j = nn - 1; j >= 0; j--)
            sumn = sumn * x + cofs.$()[j];
        for (j = nn + dd - 2; j >= nn; j--)
            sumd = sumd * x + cofs.$()[j];
        return sumn / (1.0 + x * sumd);
    }
}
