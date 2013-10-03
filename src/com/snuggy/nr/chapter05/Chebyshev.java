
package com.snuggy.nr.chapter05;

import static com.snuggy.nr.refs.Refs.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.refs.*;
import com.snuggy.nr.util.*;

@Deprecated @Broken
public class Chebyshev implements Func_Doub_To_Doub, ByValue<Chebyshev> {

    // Object for Chebyshev approximation and related methods.
    private int n, m; // Number of total, and truncated, coefficients.
    private double[] c;
    private double a, b; // Approximation interval.

    // Chebyshev(Doub func(Doub), Doub aa, Doub bb, Int nn);

    // Constructor. Approximate the function func in the interval [aa,bb]
    // with nn terms.
    public Chebyshev(final double[] cc, final double aa, final double bb) {
        n = (cc.length);
        m = (n);
        c = doub_arr(cc);
        a = (aa);
        b = (bb);
    }
    
    private Chebyshev() {
    }

    @Override
    public Chebyshev copyOut() {
        Chebyshev r = new Chebyshev();
	    r.n = n;
	    r.m = n;
	    System.arraycopy(c, 0, r.c, 0, c.length);
	    r.a = a;
	    r.b = b;
	    return r;
    }

    @Override
    public void copyIn(Chebyshev t) {
	    n = t.n;
	    m = t.n;
	    System.arraycopy(t.c, 0, c, 0, t.c.length);
	    a = t.a;
	    b = t.b;
    }

    // Constructor from previously computed coefficients.
    public int setm(double thresh) {
        while (m > 1 && abs(c[m - 1]) < thresh)
            m--;
        return m;
    }

    public double[] c() {
        return c;
    }

    // Set m, the number of coefficients after truncating to an error level
    // thresh, and return the value set.
    // Doub eval(Doub x, Int m);

    public double eval(final double x) throws NRException {
        return eval(x, m);
    }

    // Return a value for the Chebyshev fit, either using the stored m or else
    // overriding it.

    // Chebyshev derivative(); See 5.9.
    // Chebyshev integral();
    // VecDoub polycofs(Int m); See 5.10.

    public double[] polycofs() {
        return polycofs(m);
    }

    // Chebyshev(VecDoub &pc); // See 5.11.

    public Chebyshev(Func_Doub_To_Doub func, final double aa, final double bb) throws NRException {
        this(func, aa, bb, 50);
    }

    public Chebyshev(Func_Doub_To_Doub func, final double aa, final double bb, final int nn) throws NRException {
        // Chebyshev fit: Given a function func, lower and upper limits of the
        // interval [a,b], compute and save nn coefficients of the Chebyshev
        // approximation such that func.x/ OE
        // Pnn-1
        // kD0 ckTk.y/ 
        // c0=2, where y and x are related by (5.8.10). This routine is
        // intended to be called with moderately large n (e.g., 30 or 50),
        // the array of c’s subsequently to be truncated at the smaller value
        // m such that cm and subsequent elements are negligible.
        n = (nn);
        m = (nn);
        c = doub_arr(n);
        a = (aa);
        b = (bb);
        final double pi = 3.141592653589793;
        int k, j;
        double fac, bpa, bma, y, sum;
        double[] f = doub_arr(n);
        bma = 0.5 * (b - a);
        bpa = 0.5 * (b + a);
        for (k = 0; k < n; k++) { // We evaluate the function at the n points
                                  // required
            y = cos(pi * (k + 0.5) / n); // by (5.8.7).
            f[k] = func.eval(y * bma + bpa);
        }
        fac = 2.0 / n;
        for (j = 0; j < n; j++) { // Now evaluate (5.8.7).
            sum = 0.0;
            for (k = 0; k < n; k++)
                sum += f[k] * cos(pi * j * (k + 0.5) / n);
            c[j] = fac * sum;
        }
    }

    public double eval(final double x, final int m) throws NRException {
        // Chebyshev evaluation: The Chebyshev polynomial
        // Pm-1
        // kD0 ckTk.y/  c0=2 is evaluated at a point y D
        // OEx  .bCa/=2 =OE.b  a/=2 , and the result is returned as the
        // function value.
        double d = 0.0, dd = 0.0, sv, y, y2;
        int j;
        if ((x - a) * (x - b) > 0.0)
            throw new NRException("x not in range in Chebyshev::eval");
        y2 = 2.0 * (y = (2.0 * x - a - b) / (b - a)); // Change of variable.
        for (j = m - 1; j > 0; j--) { // Clenshaw’s recurrence.
            sv = d;
            d = y2 * d - dd + c[j];
            dd = sv;
        }
        return y * d - dd + 0.5 * c[0]; // Last step is different.
    }

    public $$<Chebyshev> derivative() {
        // Return a new Chebyshev object that approximates the derivative
        // of the existing function over the same range [a,b].
        int j;
        double con;
        double[] cder = doub_arr(n);
        cder[n - 1] = 0.0; // n-1 and n-2 are special cases.
        cder[n - 2] = 2 * (n - 1) * c[n - 1];
        for (j = n - 2; j > 0; j--)
            // Equation (5.9.2).
            cder[j - 1] = cder[j + 1] + 2 * j * c[j];
        con = 2.0 / (b - a);
        for (j = 0; j < n; j++)
            cder[j] *= con; // Normalize to the interval b-a.
        return $$(new Chebyshev(cder, a, b));
    }

    public $$<Chebyshev> integral() {
        // Return a new Chebyshev object that approximates the indefinite
        // integral of the existing function over the same range [a,b]. The
        // constant of integration is set so that the integral vanishes at a.
        int j;
        double sum = 0.0, fac = 1.0, con;
        double[] cint = doub_arr(n);
        con = 0.25 * (b - a); // Factor that normalizes to the interval b-a.
        for (j = 1; j < n - 1; j++) {
            cint[j] = con * (c[j - 1] - c[j + 1]) / j; // Equation (5.9.1).
            sum += fac * cint[j]; // Accumulates the constant of integration.
            fac = -fac; // Will equal ?1.
        }
        cint[n - 1] = con * c[n - 2] / (n - 1); // Special case of (5.9.1) for
                                                // n-1.
        sum += fac * cint[n - 1];
        cint[0] = 2.0 * sum; // Set the constant of integration.
        return $$(new Chebyshev(cint, a, b));
    }

    public double[] polycofs(final int m) {
        // Polynomial coefficients from a Chebyshev fit. Given a coefficient
        // array c[0..n-1], this routine returns a coefficient array d[0..n-1]
        // such that
        // Pn-1
        // kD0 dkyk D
        // Pn-1
        // kD0 ckTk.y/  c0=2. The method is Clenshaw’s recurrence (5.8.11),
        // but now applied algebraically rather than arithmetically.
        int k, j;
        double sv;
        double[] d = doub_arr(m), dd = doub_arr(m);
        for (j = 0; j < m; j++)
            d[j] = dd[j] = 0.0;
        d[0] = c[m - 1];
        for (j = m - 2; j > 0; j--) {
            for (k = m - j; k > 0; k--) {
                sv = d[k];
                d[k] = 2.0 * d[k - 1] - dd[k];
                dd[k] = sv;
            }
            sv = d[0];
            d[0] = -dd[0] + c[j];
            dd[0] = sv;
        }
        for (j = m - 1; j > 0; j--)
            d[j] = d[j - 1] - dd[j];
        d[0] = -dd[0] + 0.5 * c[0];
        return d;
    }

    public Chebyshev(final double[] d) {
        // Inverse of routine polycofs in Chebyshev: Given an array of
        // polynomial coefficients d[0..n-1], construct an equivalent
        // Chebyshev object.
        n = (d.length);
        m = (n);
        c = doub_arr(n);
        a = (-1.);
        b = (1.);
        c[n - 1] = d[n - 1];
        c[n - 2] = 2.0 * d[n - 2];
        for (int j = n - 3; j >= 0; j--) {
            c[j] = 2.0 * d[j] + c[j + 2];
            for (int i = j + 1; i < n - 2; i++) {
                c[i] = (c[i] + c[i + 2]) / 2;
            }
            c[n - 2] /= 2;
            c[n - 1] /= 2;
        }
    }
}
