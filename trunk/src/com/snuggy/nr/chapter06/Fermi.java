
package com.snuggy.nr.chapter06;

import static java.lang.Math.*;

import com.snuggy.nr.chapter04.*;
import com.snuggy.nr.util.*;

public class Fermi implements Func_Doub_To_Doub, Func_Doub_Doub_To_Doub {

    private double kk, etaa, thetaa;

    // Doub operator() (final double t);
    // Doub operator() (final double x, final double del);
    // Doub val(final double k, final double eta, final double theta);

    public double eval(final double t) {
        // Integrand for trapezoidal quadrature of generalized Fermi-Dirac
        // integral with transformation
        // x D exp.t  et /.
        double x;
        x = exp(t - exp(-t));
        return x * (1.0 + exp(-t)) * pow(x, kk) * sqrt(1.0 + thetaa * 0.5 * x) / (exp(x - etaa) + 1.0);
    }

    public double eval(final double x, final double del) {
        // Integrand for DE rule quadrature of generalized Fermi-Dirac integral.
        if (x < 1.0)
            return pow(del, kk) * sqrt(1.0 + thetaa * 0.5 * x) / (exp(x - etaa) + 1.0);
        else
            return pow(x, kk) * sqrt(1.0 + thetaa * 0.5 * x) / (exp(x - etaa) + 1.0);
    }

    public double val(final double k, final double eta, final double theta) throws NRException {
        // Computes the generalized Fermi-Dirac integral Fk.
        // ; /, where k > 1 and 
        // 0. The accuracy is approximately the square of the parameter EPS.
        // NMAX limits the total number of quadrature steps.
        final double EPS = 3.0e-9;
        final int NMAX = 11;
        double a, aa, b, bb, hmax, olds = 0.0, sum;
        kk = k; // Load the arguments into the member variables
        etaa = eta; // for use in the function evaluations.
        thetaa = theta;
        if (eta <= 15.0) {
            a = -4.5; // Set limits for x D exp.t  et / mapping.
            b = 5.0;
            Trapzd<Fermi> s = new Trapzd<Fermi>(this, a, b);
            for (int i = 1; i <= NMAX; i++) {
                sum = s.next();
                if (i > 3) // Test for convergence.
                    if (abs(sum - olds) <= EPS * abs(olds))
                        return sum;
                olds = sum; // Save value for next convergence test.
            }
        } else {
            a = 0.0; // Set limits for DE rule.
            b = eta;
            aa = eta;
            bb = eta + 60.0;
            hmax = 4.3; // Big enough to handle negative k or large .
            DErule<Fermi> s = new DErule<Fermi>(this, a, b, hmax);
            DErule<Fermi> ss = new DErule<Fermi>(this, aa, bb, hmax);
            for (int i = 1; i <= NMAX; i++) {
                sum = s.next() + ss.next();
                if (i > 3)
                    if (abs(sum - olds) <= EPS * abs(olds))
                        return sum;
                olds = sum;
            }
        }
        throw new NRException("no convergence in fermi");
        // return 0.0;
    }

}
