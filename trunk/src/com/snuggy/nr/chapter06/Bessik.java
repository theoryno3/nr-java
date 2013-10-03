
package com.snuggy.nr.chapter06;

import static com.snuggy.nr.refs.Refs.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.refs.*;

public class Bessik {

    private static final double[] i0p = doub_arr(14), i0q = doub_arr(5), i0pp = doub_arr(5), i0qq = doub_arr(6);
    private static final double[] i1p = doub_arr(14), i1q = doub_arr(5), i1pp = doub_arr(5), i1qq = doub_arr(6);
    private static final double[] k0pi = doub_arr(5), k0qi = doub_arr(3), k0p = doub_arr(5), k0q = doub_arr(3),
            k0pp = doub_arr(8), k0qq = doub_arr(8);
    private static final double[] k1pi = doub_arr(5), k1qi = doub_arr(3), k1p = doub_arr(5), k1q = doub_arr(3),
            k1pp = doub_arr(8), k1qq = doub_arr(8);

    private double y, z, ax, term;

    public double i0(final double x) {
        // Returns the modified Bessel function I0.x/ for any real x.
        if ((ax = abs(x)) < 15.0) { // Rational approximation.
            y = x * x;
            return poly(i0p, 0, 13, y) / poly(i0q, 0, 4, 225. - y);
        } else { // Rational approximation with ex=
            // p
            // x factored out.
            z = 1.0 - 15.0 / ax;
            return exp(ax) * poly(i0pp, 0, 4, z) / (poly(i0qq, 0, 5, z) * sqrt(ax));
        }
    }

    public double i1(final double x) {
        // Returns the modified Bessel function I1.x/ for any real x.
        if ((ax = abs(x)) < 15.0) { // Rational approximation.
            y = x * x;
            return x * poly(i1p, 0, 13, y) / poly(i1q, 0, 4, 225. - y);
        } else { // Rational approximation with ex=
            // p
            // x factored out.
            z = 1.0 - 15.0 / ax;
            double ans = exp(ax) * poly(i1pp, 0, 4, z) / (poly(i1qq, 0, 5, z) * sqrt(ax));
            return x > 0.0 ? ans : -ans;
        }
    }

    public double k0(final double x) {
        // Returns the modified Bessel function K0.x/ for positive real x.
        if (x <= 1.0) { // Use two rational approximations.
            z = x * x;
            term = poly(k0pi, 0, 4, z) * log(x) / poly(k0qi, 0, 2, 1. - z);
            return poly(k0p, 0, 4, z) / poly(k0q, 0, 2, 1. - z) - term;
        } else { // Rational approximation with ex=
            // p
            // x factored
            z = 1.0 / x; // out.
            return exp(-x) * poly(k0pp, 0, 7, z) / (poly(k0qq, 0, 7, z) * sqrt(x));
        }
    }

    public double k1(final double x) {
        // Returns the modified Bessel function K1.x/ for positive real x.
        if (x <= 1.0) { // Use two rational approximations.
            z = x * x;
            term = poly(k1pi, 0, 4, z) * log(x) / poly(k1qi, 0, 2, 1. - z);
            return x * (poly(k1p, 0, 4, z) / poly(k1q, 0, 2, 1. - z) + term) + 1. / x;
        } else { // Rational approximation with ex=
            // p
            // x factored
            z = 1.0 / x; // out.
            return exp(-x) * poly(k1pp, 0, 7, z) / (poly(k1qq, 0, 7, z) * sqrt(x));
        }
    }

    // Doub in(final int n, final double x);
    // Returns the modified Bessel function In.x/ for any real x and n 0.
    // Doub kn(final int n, final double x);
    // Returns the modified Bessel function Kn.x/ for positive x and n 0.

    public double poly(final double[] cof_arr, final int cof_off, final int n, final double x) {
        // Common code: Evaluate a polynomial.
        double ans = cof_arr[cof_off + n];
        for (int i = n - 1; i >= 0; i--)
            ans = ans * x + cof_arr[cof_off + i];
        return ans;
    }

    public double kn(final int n, final double x) {
        // Returns the modified Bessel function Kn.x/ for positive x and n 0.
        int j;
        double bk, bkm, bkp, tox;
        if (n == 0)
            return k0(x);
        if (n == 1)
            return k1(x);
        tox = 2.0 / x;
        bkm = k0(x); // Upward recurrence for all x...
        bk = k1(x);
        for (j = 1; j < n; j++) { // ...and here it is.
            bkp = bkm + j * tox * bk;
            bkm = bk;
            bk = bkp;
        }
        return bk;
    }

    public double in(final int n, final double x) {
        // Returns the modified Bessel function In.x/ for any real x and n 0.
        final double ACC = 200.0; // ACC determines accuracy.
        final int IEXP = Double.MAX_EXPONENT / 2; // numeric_limits<Doub>::max_exponent/2;
        int j;
        $int k = $(0);
        double bi, bim, bip, /* dum, */tox, ans;
        if (n == 0)
            return i0(x);
        if (n == 1)
            return i1(x);
        if (x * x <= 8.0 * Double.MIN_VALUE) // numeric_limits<Doub>::min())
            return 0.0;
        else {
            tox = 2.0 / abs(x);
            bip = ans = 0.0;
            bi = 1.0;
            for (j = 2 * (n + Int(sqrt(ACC * n))); j > 0; j--) { // Downward
                                                                 // recurrence.
                bim = bip + j * tox * bi;
                bip = bi;
                bi = bim;
                /* dum = */frexp(bi, k);
                if (k.$() > IEXP) { // Renormalize to prevent overflows.
                    ans = ldexp(ans, -IEXP);
                    bi = ldexp(bi, -IEXP);
                    bip = ldexp(bip, -IEXP);
                }
                if (j == n)
                    ans = bip;
            }
            ans *= i0(x) / bi; // Normalize with bessi0.
            return (x < 0.0) && ((n & 1) != 0) ? -ans : ans;
        }
    }
}
