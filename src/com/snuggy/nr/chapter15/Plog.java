package com.snuggy.nr.chapter15;

import static com.snuggy.nr.chapter06.Static.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Plog {

    // Functor that calculates logP of a State.
    private final double[] dat; // Bind to data vector.
    private int ndat;
    private final double[] stau, slogtau;

    public Plog(final double[] data) {
        // Constructor. Digest the data vector for subsequent fast calculation
        // of logP. The data are assumed to be sorted in ascending order.
        dat = (data);
        ndat = (data.length);
        stau = doub_vec(ndat);
        slogtau = doub_vec(ndat);
        int i;
        stau[0] = slogtau[0] = 0.;
        for (i = 1; i < ndat; i++) {
            stau[i] = dat[i] - dat[0]; // Equal to sum of intervals.
            slogtau[i] = slogtau[i - 1] + log(dat[i] - dat[i - 1]);
        }
    }

    public double func(final State s) throws NRException {
        // Return logP of s, and also set s.plog.
        int i, ilo, ihi, n1, n2;
        double st1, st2, stl1, stl2, ans;
        ilo = 0;
        ihi = ndat - 1;
        while (ihi - ilo > 1) { // Bisection to find where is tc in the data.
            i = (ihi + ilo) >> 1;
            if (s.tc() > dat[i])
                ilo = i;
            else
                ihi = i;
        }
        n1 = ihi;
        n2 = ndat - 1 - ihi;
        st1 = stau[ihi];
        st2 = stau[ndat - 1] - st1;
        stl1 = slogtau[ihi];
        stl2 = slogtau[ndat - 1] - stl1;
        ans = n1 * (s.k1() * log(s.lam1()) - factln(s.k1() - 1)) + (s.k1() - 1) * stl1 - s.lam1() * st1;
        ans += n2 * (s.k2() * log(s.lam2()) - factln(s.k2() - 1)) + (s.k2() - 1) * stl2 - s.lam2() * st2;
        s.set_plog(ans);
        return ans;
    }

}
