
package com.snuggy.nr.chapter15;

import static java.lang.Math.*;

import com.snuggy.nr.chapter07.*;
import com.snuggy.nr.refs.*;

public class Proposal {

    // Functor implementing the proposal distribution.
    private Normaldev gau;
    private double logstep;

    public Proposal(final int ranseed, final double lstep) {
        gau = new Normaldev(0., 1., ranseed);
        logstep = (lstep);
    }
    
    public Normaldev gau() {
        return gau;
    }
    
    //private static long hit_count = 0;

    public void func(final State s1, final State s2, final $double qratio) {
        // Given state s1, set state s2 to a proposed candidate. Also set
        // qratio to q.s1js2/=q.s2js1/.
        //hit_count++;
        //if (hit_count == 218)
	    //    System.out.println("Proposal.hit_count is " + hit_count);
        double r = gau.doub();
        if (r < 0.9) { // Lognormal steps holding the k’s constant.
            s2.set_lam1(s1.lam1() * exp(logstep * gau.dev()));
            s2.set_lam2(s1.lam2() * exp(logstep * gau.dev()));
            s2.set_tc(s1.tc() * exp(logstep * gau.dev()));
            s2.set_k1(s1.k1());
            s2.set_k2(s1.k2());
            qratio.$((s2.lam1() / s1.lam1()) * (s2.lam2() / s1.lam2()) * (s2.tc() / s1.tc()));
            // Factors for lognormal steps.
        } else { // Steps that change k1 and/or k2.
            r = gau.doub();
            if (s1.k1() > 1) {
                if (r < 0.5)
                    s2.set_k1(s1.k1());
                else if (r < 0.75)
                    s2.set_k1(s1.k1() + 1);
                else
                    s2.set_k1(s1.k1() - 1);
            } else { // k1 D 1 requires special treatment.
                if (r < 0.75)
                    s2.set_k1(s1.k1());
                else
                    s2.set_k1(s1.k1() + 1);
            }
            s2.set_lam1(s2.k1() * s1.lam1() / s1.k1());
            r = gau.doub(); // Now all the same for k2.
            if (s1.k2() > 1) {
                if (r < 0.5)
                    s2.set_k2(s1.k2());
                else if (r < 0.75)
                    s2.set_k2(s1.k2() + 1);
                else
                    s2.set_k2(s1.k2() - 1);
            } else {
                if (r < 0.75)
                    s2.set_k2(s1.k2());
                else
                    s2.set_k2(s1.k2() + 1);
            }
            s2.set_lam2(s2.k2() * s1.lam2() / s1.k2());
            s2.set_tc(s1.tc());
            qratio.$(1.);
        }
    }

}
