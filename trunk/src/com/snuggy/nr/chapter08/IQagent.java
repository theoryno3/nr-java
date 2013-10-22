
package com.snuggy.nr.chapter08;

import static com.snuggy.nr.chapter08.Static.*;
import static com.snuggy.nr.refs.Refs.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.refs.*;
import com.snuggy.nr.util.*;

public class IQagent {

    // Object for estimating arbitrary quantile values from a continuing
    // stream of data values.
    private static final int nbuf = 1000; // Batch size. You may 10 if you
                                          // expect >
    private int nq, nt, nd; // 106 data values.
    private double q0, qm;
    private final double[] pval, dbuf;
    private final $$double1d qile;

    public IQagent() throws NRException {
        nq = (251);
        nt = (0);
        nd = (0);
        pval = doub_vec(nq);
        dbuf = doub_vec(nbuf);
        qile = $$(doub_vec(nq, 0.));
        q0 = (1.e99);
        qm = (-1.e99);
        // Constructor. No arguments.
        for (int j = 85; j <= 165; j++)
            pval[j] = (j - 75.) / 100.;
        // Set general purpose array of p-values ranging from 106 to 1106.
        // You can change this if you want:
        for (int j = 84; j >= 0; j--) {
            pval[j] = 0.87191909 * pval[j + 1];
            pval[250 - j] = 1. - pval[j];
        }
    }

    public void add(final double datum) throws NRException {
        // Assimilate a new value from the stream.
        dbuf[nd++] = datum;
        if (datum < q0) {
            q0 = datum;
        }
        if (datum > qm) {
            qm = datum;
        }
        if (nd == nbuf)
            update(); // Time for a batch update.
    }

    public void update() throws NRException {
        // Batch update, as shown in Figure 8.5.1. This function is called by
        // add or report and should not be called directly by the user.
        int jd = 0, jq = 1, iq;
        double target, told = 0., tnew = 0., qold, qnew;
        double[] newqile = doub_vec(nq); // Will be new quantiles after update.
        sort(dbuf, nd);
        qold = qnew = qile.$()[0] = newqile[0] = q0; // Set lowest and highest
                                                     // to min
        // and max values seen so far, and set compatible p-values.
        qile.$()[nq - 1] = newqile[nq - 1] = qm;
        pval[0] = min(0.5 / (nt + nd), 0.5 * pval[1]);
        pval[nq - 1] = max(1. - 0.5 / (nt + nd), 0.5 * (1. + pval[nq - 2]));
        for (iq = 1; iq < nq - 1; iq++) { // Main loop over target p-values for
                                          // inter
            target = (nt + nd) * pval[iq]; // polation.
            if (tnew < target)
                for (;;) {
                    // Here’s the guts: We locate a succession of
                    // abscissa-ordinate
                    // pairs (qnew,tnew) that are the discontinuities of value
                    // or slope
                    // in Figure 8.5.1(c), breaking to perform an interpolation
                    // as we cross
                    // each target.
                    if (jq < nq && (jd >= nd || qile.$()[jq] < dbuf[jd])) {
                        // Found slope discontinuity from old CDF.
                        qnew = qile.$()[jq];
                        tnew = jd + nt * pval[jq++];
                        if (tnew >= target)
                            break;
                    } else { // Found value discontinuity from batch data
                        qnew = dbuf[jd]; // CDF.
                        tnew = told;
                        if (qile.$()[jq] > qile.$()[jq - 1])
                            tnew += nt * (pval[jq] - pval[jq - 1]) * (qnew - qold) / (qile.$()[jq] - qile.$()[jq - 1]);
                        jd++;
                        if (tnew >= target)
                            break;
                        told = tnew++;
                        qold = qnew;
                        if (tnew >= target)
                            break;
                    }
                    told = tnew;
                    qold = qnew;
                } // Break to here and perform the new interpolation.
            if (tnew == told)
                newqile[iq] = 0.5 * (qold + qnew);
            else
                newqile[iq] = qold + (qnew - qold) * (target - told) / (tnew - told);
            told = tnew;
            qold = qnew;
        }
        // qile = newqile;
        $$(qile, newqile);
        nt += nd;
        nd = 0;
    }

    public double report(final double p) throws NRException {
        // Return estimated p-quantile for the data seen so far. (E.g., p D 0:5
        // for median.)
        double q;
        if (nd > 0)
            update(); // You may want to remove this line. See text.
        int jl = 0, jh = nq - 1, j;
        while (jh - jl > 1) { // Locate place in table by bisection.
            j = (jh + jl) >> 1;
            if (p > pval[j])
                jl = j;
            else
                jh = j;
        }
        j = jl; // Interpolate.
        q = qile.$()[j] + (qile.$()[j + 1] - qile.$()[j]) * (p - pval[j]) / (pval[j + 1] - pval[j]);
        return MAX(qile.$()[0], MIN(qile.$()[nq - 1], q));
    }

}
