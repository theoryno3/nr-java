
package com.snuggy.nr.chapter06;

import static com.snuggy.nr.chapter06.Static.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Gamma extends Gauleg18 {
    // Object for incomplete gamma function. Gauleg18 provides coefficients
    // for Gauss-Legendre quadrature.
    private static final int ASWITCH = 100; // When to switch to quadrature
                                            // method.
    // private static final double EPS; // See end of struct for
    // initializations.
    // private static final double FPMIN;
    private static final double EPS = EPS(); // numeric_limits<Doub>::epsilon();
    private static final double FPMIN = Double.MIN_VALUE / EPS; // numeric_limits<Doub>::min()/EPS;

    private double gln;

    public double gammp(final double a, final double x) throws NRException {
        // Returns the incomplete gamma function P.a;x/.
        if (x < 0.0 || a <= 0.0)
            throw new NRException("bad args in gammp");
        if (x == 0.0)
            return 0.0;
        else if ((int) a >= ASWITCH)
            return gammpapprox(a, x, 1); // Quadrature.
        else if (x < a + 1.0)
            return gser(a, x); // Use the series representation.
        else
            return 1.0 - gcf(a, x); // Use the continued fraction
                                    // representation.
    }

    public double gammq(final double a, final double x) throws NRException {
        // Returns the incomplete gamma function Q.a; x/  1  P.a;x/.
        if (x < 0.0 || a <= 0.0)
            throw new NRException("bad args in gammq");
        if (x == 0.0)
            return 1.0;
        else if ((int) a >= ASWITCH)
            return gammpapprox(a, x, 0); // Quadrature.
        else if (x < a + 1.0)
            return 1.0 - gser(a, x); // Use the series representation.
        else
            return gcf(a, x); // Use the continued fraction representation.
    }

    public double gser(final double a, final double x) throws NRException {
        // Returns the incomplete gamma function P.a;x/ evaluated by its
        // series representation.
        // Also sets ln
        // .a/ as gln. User should not call directly.
        double sum, del, ap;
        gln = gammln(a);
        ap = a;
        del = sum = 1.0 / a;
        for (;;) {
            ++ap;
            del *= x / ap;
            sum += del;
            if (abs(del) < abs(sum) * EPS) {
                return sum * exp(-x + a * log(x) - gln);
            }
        }
    }

    public double gcf(final double a, final double x) throws NRException {
        // Returns the incomplete gamma function Q.a; x/ evaluated by its
        // continued fraction representation. Also sets ln
        // .a/ as gln. User should not call directly.
        int i;
        double an, b, c, d, del, h;
        gln = gammln(a);
        b = x + 1.0 - a; // Set up for evaluating continued fraction
        // by modified Lentz’s method (5.2) with b0 D 0.
        c = 1.0 / FPMIN;
        d = 1.0 / b;
        h = d;
        for (i = 1;; i++) { // Iterate to convergence.
            an = -i * (i - a);
            b += 2.0;
            d = an * d + b;
            if (abs(d) < FPMIN)
                d = FPMIN;
            c = b + an / c;
            if (abs(c) < FPMIN)
                c = FPMIN;
            d = 1.0 / d;
            del = d * c;
            h *= del;
            if (abs(del - 1.0) <= EPS)
                break;
        }
        return exp(-x + a * log(x) - gln) * h; // Put factors in front.
    }

    public double gammpapprox(final double a, final double x, final int psig) throws NRException {
        // Incomplete gamma by quadrature. Returns P.a;x/ or Q.a; x/, when
        // psig is 1 or 0, respectively. User should not call directly.
        int j;
        double xu, t, sum, ans;
        double a1 = a - 1.0, lna1 = log(a1), sqrta1 = sqrt(a1);
        gln = gammln(a);
        // Set how far to integrate into the tail:
        if (x > a1)
            xu = MAX(a1 + 11.5 * sqrta1, x + 6.0 * sqrta1);
        else
            xu = MAX(0., MIN(a1 - 7.5 * sqrta1, x - 5.0 * sqrta1));
        sum = 0;
        for (j = 0; j < ngau; j++) { // Gauss-Legendre.
            t = x + (xu - x) * y[j];
            sum += w[j] * exp(-(t - a1) + a1 * (log(t) - lna1));
        }
        ans = sum * (xu - x) * exp(a1 * (lna1 - 1.) - gln);
        return (psig != 0 ? (ans > 0.0 ? 1.0 - ans : -ans) : (ans >= 0.0 ? ans : 1.0 + ans));
    }

    // Doub invgammp(Doub p, Doub a);
    // Inverse function on x of P.a;x/. See 6.2.1.

    public double invgammp(final double p, final double a) throws NRException {
        // Returns x such that P.a;x/ D p for an argument p between 0 and 1.
        int j;
        double x, err, t, u, pp, lna1 = 0.0, afac = 0.0, a1 = a - 1;
        final double EPS = 1.e-8; // Accuracy is the square of EPS.
        gln = gammln(a);
        if (a <= 0.)
            throw new NRException("a must be pos in invgammap");
        if (p >= 1.)
            return MAX(100., a + 100. * sqrt(a));
        if (p <= 0.)
            return 0.0;
        if (a > 1.) { // Initial guess based on reference [1].
            lna1 = log(a1);
            afac = exp(a1 * (lna1 - 1.) - gln);
            pp = (p < 0.5) ? p : 1. - p;
            t = sqrt(-2. * log(pp));
            x = (2.30753 + t * 0.27061) / (1. + t * (0.99229 + t * 0.04481)) - t;
            if (p < 0.5)
                x = -x;
            x = MAX(1.e-3, a * pow(1. - 1. / (9. * a) - x / (3. * sqrt(a)), 3));
        } else { // Initial guess based on equations (6.2.8)
            t = 1.0 - a * (0.253 + a * 0.12); // and (6.2.9).
            if (p < t)
                x = pow(p / t, 1. / a);
            else
                x = 1. - log(1. - (p - t) / (1. - t));
        }
        for (j = 0; j < 12; j++) {
            if (x <= 0.0)
                return 0.0; // x too small to compute accurately.
            err = gammp(a, x) - p;
            if (a > 1.)
                t = afac * exp(-(x - a1) + a1 * (log(x) - lna1));
            else
                t = exp(-x + a1 * log(x) - gln);
            u = err / t;
            x -= (t = u / (1. - 0.5 * MIN(1., u * ((a - 1.) / x - 1)))); // Halley’s
                                                                         // method.
            if (x <= 0.)
                x = 0.5 * (x + t); // Halve old value if x tries to go negative.
            if (abs(t) < EPS * x)
                break;
        }
        return x;
    }
}
