
package com.snuggy.nr.chapter07;

import static com.snuggy.nr.refs.Refs.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.refs.*;
import com.snuggy.nr.util.*;

public class MCintegrate {

    // Object for Monte Carlo integration of one or more functions in an
    // ndim-dimensional region.
    private int ndim, nfun, n; // Number of dimensions, functions, and points
                               // sampled.
    public final double[] ff, fferr; // Answers: The integrals and their standard
                                // errors.
    private final $$double1d xlo, xhi;
    private final double[] x, sf, sferr;
    private final double[] xx, fn;
    private double vol; // Volume of the box V .
    
    private Func_DoubArr_To_DoubArr funcsp; // Pointers to the user-supplied functions.
    private Func_DoubArr_To_DoubArr xmapp;
    private Func_DoubArr_To_Bool inregionp;
    private Ran ran; // Random number generator.
    
    // MCintegrate(const VecDoub &xlow, const VecDoub &xhigh,
    // VecDoub funcs(const VecDoub &), Bool inregion(const VecDoub &),
    // VecDoub xmap(const VecDoub &), Int ranseed);
    // Constructor. The arguments are in the order described in the itemized
    // list above.
    // void step(Int nstep);
    // Sample an additional nstep points, accumulating the various sums.
    // void calcanswers();

    public MCintegrate(final double[] xlow, final double[] xhigh, 
                        final Func_DoubArr_To_DoubArr funcs, final Func_DoubArr_To_Bool inregion,
                        final Func_DoubArr_To_DoubArr xmap, final int ranseed) throws NRException {
        ndim = (xlow.length);
        n = (0);
        xlo = $$(xlow);
        xhi = $$(xhigh);
        x = (doub_vec(ndim));
        xx = (new double[ndim]);
        funcsp = (funcs);
        xmapp = (xmap);
        inregionp = (inregion);
        vol = (1.);
        ran = new Ran(ranseed);
        if (xmapp != null)
            nfun = funcs.eval(xmapp.eval(xlo.$())).length;
        else
            nfun = funcs.eval(xlo.$()).length;
        ff = (doub_vec(nfun));
        fferr = (doub_vec(nfun));
        fn = (new double[nfun]);
        sf = (doub_vec(nfun, 0.));
        sferr = (doub_vec(nfun, 0.));
        for (int j = 0; j < ndim; j++)
            vol *= abs(xhi.$()[j] - xlo.$()[j]);
    }

    public void step(final int nstep) throws NRException {
        int i, j;
        for (i = 0; i < nstep; i++) {
            for (j = 0; j < ndim; j++)
                x[j] = xlo.$()[j] + (xhi.$()[j] - xlo.$()[j]) * ran.doub();
            if (xmapp != null)
                $$(xx, xmapp.eval(x));
            else
                $$(xx, x);
            if (inregionp.eval(xx)) {
                $$(fn, funcsp.eval(xx));
                for (j = 0; j < nfun; j++) {
                    sf[j] += fn[j];
                    sferr[j] += SQR(fn[j]);
                }
            }
        }
        n += nstep;
    }

    public void calcanswers() {
        for (int j = 0; j < nfun; j++) {
            ff[j] = vol * sf[j] / n;
            fferr[j] = vol * sqrt((sferr[j] / n - SQR(sf[j] / n)) / n);
        }
    }
}
