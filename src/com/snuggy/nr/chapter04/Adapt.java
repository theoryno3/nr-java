package com.snuggy.nr.chapter04;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Adapt {

    // Adaptive quadrature for the integral of the function func from a to b.
    // Integration is performed by a Gauss-Lobatto method with a
    // Kronrod extension.
    // private static final double alpha,beta,x1,x2,x3; // Abscissas for
    // Gauss-Lobatto-Kronrod
    // private static final double[] x = doub_arr(12);

    private double TOL, toler;
    private boolean terminate ; // quadrature.
    @SuppressWarnings("unused")
    private boolean out_of_tolerance; 

    private final double alpha = sqrt(2.0 / 3.0);
    private final double beta = 1.0 / sqrt(5.0);
    private final double x1 = 0.942882415695480;
    private final double x2 = 0.641853342345781;
    private final double x3 = 0.236383199662150;
    private final double x[] = { 0, -x1, -alpha, -x2, -beta, -x3, 0.0, x3, beta, x2, alpha, x1 };

    // Adapt(Doub tol);
    // template <class T>
    // Doub integrate(T &func, final double a, final double b);
    // template <class T>
    // Doub adaptlob(T &func, final double a, final double b, final double fa,
    // final double fb, final double is);
    private final double EPS = EPS(); // numeric_limits<Doub>::epsilon();

    public Adapt(final double tol) {
        // Constructor is invoked with desired tolerance tol. The smallest
        // allowed value of tol is 10*EPS, where EPS is the machine precision.
        // If tol is input as less than this (e.g. tol D 0), then tol is set
        // to 10*EPS.
        TOL = (tol);
        terminate = (true);
        out_of_tolerance = (false);
        if (TOL < 10.0 * EPS)
            TOL = 10.0 * EPS;
    }

    public <T extends Func_Doub_To_Doub> double integrate(final T func, final double a, final double b) throws NRException {
        double m, h, fa, fb, i1, i2, is, erri1, erri2, r;
        final double[] y = doub_vec(13);
        m = 0.5 * (a + b);
        h = 0.5 * (b - a);
        fa = y[0] = func.eval(a);
        fb = y[12] = func.eval(b);
        for (int i = 1; i < 12; i++)
            y[i] = func.eval(m + x[i] * h);
        i2 = (h / 6.0) * (y[0] + y[12] + 5.0 * (y[4] + y[8])); // 4-point
                                                               // Gauss-Lobatto
                                                               // formula.
        i1 = (h / 1470.0) * (77.0 * (y[0] + y[12]) + 432.0 * (y[2] + y[10]) + 625.0 * (y[4] + y[8]) + 672.0 * y[6]); // 7-point
                                                                                                                     // Kronrod
                                                                                                                     // extension.
        is = h
                * (0.0158271919734802 * (y[0] + y[12]) + 0.0942738402188500 * (y[1] + y[11]) + 0.155071987336585
                        * (y[2] + y[10]) + 0.188821573960182 * (y[3] + y[9]) + 0.199773405226859 * (y[4] + y[8])
                        + 0.224926465333340 * (y[5] + y[7]) + 0.242611071901408 * y[6]); // 13-point
                                                                                         // Kronrod
                                                                                         // extension.
        erri1 = abs(i1 - is);
        erri2 = abs(i2 - is);
        r = (erri2 != 0.0) ? erri1 / erri2 : 1.0;
        toler = (r > 0.0 && r < 1.0) ? TOL / r : TOL; // Error of i1 will be
                                                      // suciently small
        if (is == 0.0) // that we can increase tolerance.
            is = b - a;
        is = abs(is);
        return adaptlob(func, a, b, fa, fb, is);
    }

    public <T extends Func_Doub_To_Doub> double adaptlob(final T func, final double a, final double b, final double fa,
            final double fb, final double is) throws NRException {
        // Helper function for recursion.
        double m, h, mll, ml, mr, mrr, fmll, fml, fm, fmrr, fmr, i1, i2;
        m = 0.5 * (a + b);
        h = 0.5 * (b - a);
        mll = m - alpha * h;
        ml = m - beta * h;
        mr = m + beta * h;
        mrr = m + alpha * h;
        fmll = func.eval(mll);
        fml = func.eval(ml);
        fm = func.eval(m);
        fmr = func.eval(mr);
        fmrr = func.eval(mrr);
        i2 = h / 6.0 * (fa + fb + 5.0 * (fml + fmr)); // 4-point Gauss-Lobatto
                                                      // formula.
        i1 = h / 1470.0 * (77.0 * (fa + fb) + 432.0 * (fmll + fmrr) + 625.0 * (fml + fmr) + 672.0 * fm);
        // 7-point Kronrod extension.
        if (abs(i1 - i2) <= toler * is || mll <= a || b <= mrr) {
            if ((mll <= a || b <= mrr) && terminate) {
                out_of_tolerance = true; // Interval contains no more machine
                terminate = false; // numbers.
            }
            return i1; // Terminate recursion.
        } else
            // Subdivide interval.
            return adaptlob(func, a, mll, fa, fmll, is) + adaptlob(func, mll, ml, fmll, fml, is)
                    + adaptlob(func, ml, m, fml, fm, is) + adaptlob(func, m, mr, fm, fmr, is)
                    + adaptlob(func, mr, mrr, fmr, fmrr, is) + adaptlob(func, mrr, b, fmrr, fb, is);
    }
}
