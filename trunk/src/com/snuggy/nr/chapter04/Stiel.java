package com.snuggy.nr.chapter04;

import static com.snuggy.nr.chapter04.Static.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Stiel {

    // Structure for calculating the abscissas and weights of an n-point
    // Gaussian quadrature formula using the Stieltjes procedure.
    class pp implements Func_Doub_To_Doub, Func_Doub_Doub_To_Doub {
        // Functor returning the integrand for bj in eq. (4.6.7).
        private Stiel st;

        public double eval(final double x, final double del) {
            // Returns W.x/p2.x/.
            double pval = st.p(x);
            return pval * st.wt1.eval(x, del) * pval;
            // This order lessens the chance of over ow.
        }

        public double eval(final double t) throws NRException {
            // Returns W.x/p2.x/ dx=dt .
            double x = st.fx.eval(t);
            double pval = st.p(x);
            return pval * st.wt2.eval(x) * st.fdxdt.eval(t) * pval;
            // This order lessens the chance of over ow.
        }
    }

    class ppx implements Func_Doub_To_Doub, Func_Doub_Doub_To_Doub {
        // Functor returning the integrand for aj in eq. (4.6.7).
        private Stiel st;

        public double eval(final double x, final double del) {
            // Returns W.x/ x p2.x/.
            return st.ppfunc.eval(x, del) * x;
        }

        public double eval(final double t) throws NRException {
            // Returns W.x/ x p2.x/ dx=dt .
            return st.ppfunc.eval(t) * st.fx.eval(t);
        }
    }

    private int j, n;
    // Degree of the current polynomial during the computation and desired
    // number of abscissas.
    private double aa, bb, hmax;
    // For the finite range case, limits of integration and limit of integration
    // in transformed variable for the DE rule. For in nite range case, limits
    // of integration in transformed variable (hmax not used.
    private Func_Doub_Doub_To_Doub wt1; // Pointers to user-supplied functions.
    private Func_Doub_To_Doub wt2;
    private Func_Doub_To_Doub fx;
    private Func_Doub_To_Doub fdxdt;
    private final double[] a, b;
    // Coecents of the recurrence relation for the orthogonal polynomials.
    private Quadrature s1, s2;
    // The two quadratures required in each iteration of eq. (4.6.7).
    // private Doub p(final double x);
    private pp ppfunc = new pp();
    private ppx ppxfunc = new ppx();

    // Stiel(Int nn, Doub aaa, Doub bbb, Doub hmaxx, Doub wwt1(Doub,Doub));
    // Stiel(Int nn, Doub aaa, Doub bbb, Doub wwt2(Doub), Doub ffx(Doub),
    // Doub ffdxdt(Doub));
    // Doub quad(Quadrature *s);
    // void get_weights(VecDoub_O &x, VecDoub_O &w);

    public double p(final double x) {
        // Returns the orthogonal polynomial pj .x/.
        double pval = 0.0, pj, pjm1;
        if (j == 0)
            return 1.0;
        else { // Compute pj .x/ using recurrence relation.
            pjm1 = 0.0;
            pj = 1.0;
            for (int i = 0; i < j; i++) {
                pval = (x - a[i]) * pj - b[i] * pjm1;
                pjm1 = pj;
                pj = pval;
            }
        }
        return pval;
    }

    public Stiel(final int nn, final double aaa, final double bbb, final double hmaxx, Func_Doub_Doub_To_Doub wwt1) {
        // Constructor for nite-range case. Input are nn, the number of
        // quadrature abscissas and weights desired, aaa and bbb, the lower
        // and upper limits of integration, the parameter hmax to be passed to
        // the DE rule (see 4.5), and the weight function coded as a
        // function W.x; /.
        n = (nn);
        aa = (aaa);
        bb = (bbb);
        hmax = (hmaxx);
        wt1 = (wwt1);
        a = doub_vec(nn);
        b = doub_vec(nn);
        ppfunc.st = this;
        ppxfunc.st = this;
        s1 = new DErule<pp>(ppfunc, aa, bb, hmax);
        s2 = new DErule<ppx>(ppxfunc, aa, bb, hmax);
    }

    public Stiel(final int nn, final double aaa, final double bbb, Func_Doub_To_Doub wwt2, Func_Doub_To_Doub ffx, Func_Doub_To_Doub ffdxdt) {
        // Constructor for in nite-range case. Input are nn, the number of
        // quadrature abscissas and weights desired, aaa and bbb, the lower
        // and upper limits of integration, the weight function W.x/, the
        n = (nn);
        aa = (aaa);
        bb = (bbb);
        a = doub_vec(nn);
        b = doub_vec(nn);
        wt2 = (wwt2);
        fx = (ffx);
        fdxdt = (ffdxdt);
        // mapping function x.t / and its derivative dx=dt .
        ppfunc.st = this;
        ppxfunc.st = this;
        s1 = new Trapzd<pp>(ppfunc, aa, bb);
        s2 = new Trapzd<ppx>(ppxfunc, aa, bb);
    }

    public double quad(final Quadrature s) throws NRException {
        // Carries out the quadrature.
        final double EPS = 3.0e-11, MACHEPS = EPS(); // numeric_limits<Doub>::epsilon();
        // The accuracy of the quadrature is very roughly the square of EPS.
        // This choice ensures full double precision.
        final int NMAX = 11;
        double olds = 0.0, sum;
        s.n = 0;
        for (int i = 1; i <= NMAX; i++) {
            sum = s.next();
            if (i > 3)
                // Test for convergence. Modify to test absolute error if
                // integral can
                // be zero.
                if (abs(sum - olds) <= EPS * abs(olds))
                    return sum;
            if (i == NMAX)
                if (abs(sum) <= MACHEPS && abs(olds) <= MACHEPS)
                    return 0.0;
            olds = sum;
        }
        throw new NRException("no convergence in quad");
        // return 0.0;
    }

    public void get_weights(final double[] x, final double[] w) throws NRException {
        // This function returns arrays x[0..n-1] and w[0..n-1] of length n,
        // containing the abscissas and weights of the n-point Gaussian
        // quadrature formula for the weight function W.x/.
        double amu0, c, oldc = 1.0;
        if (n != x.length)
            throw new NRException("bad array size in Stiel");
        for (int i = 0; i < n; i++) { // Compute a and b arrays via eq. (4.6.7).
            j = i; // Keep track of j, degree of current polynomial.
            c = quad(s1);
            b[i] = c / oldc; // Use b[0] to store 0 D
            // R
            // W.x/ dx.
            a[i] = quad(s2) / c;
            oldc = c;
        } // The coecients a and b of the recurrence relation
          // are available at this point if needed for other purposes.
        amu0 = b[0];
        gaucof(a, b, amu0, x, w);
    }
}
