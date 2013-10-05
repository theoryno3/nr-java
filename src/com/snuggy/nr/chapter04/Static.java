
package com.snuggy.nr.chapter04;

import static com.snuggy.nr.chapter06.Static.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.chapter03.*;
import com.snuggy.nr.chapter11.*;
import com.snuggy.nr.util.*;

public class Static {

    public static <T extends Func_Doub_To_Doub> double qtrap(final T func, final double a, final double b) throws NRException {
        return qtrap(func, a, b, 1.0e-10);
    }

    public static <T extends Func_Doub_To_Doub> double qtrap(final T func, final double a, final double b, final double eps)
            throws NRException {
        // Returns the integral of the function or functor func from a to b.
        // The constants EPS can be set to the desired fractional accuracy and
        // JMAX so that 2 to the power JMAX-1 is the maximum allowed number of
        // steps. Integration is performed by the trapezoidal rule.

        final int JMAX = 20;
        double s, olds = 0.0; // Initial value of olds is arbitrary.
        Trapzd<T> t = new Trapzd<T>(func, a, b);
        for (int j = 0; j < JMAX; j++) {
            s = t.next();
            if (j > 5) // Avoid spurious early convergence.
                if (abs(s - olds) < eps * abs(olds) || (s == 0.0 && olds == 0.0))
                    return s;
            olds = s;
        }
        throw new NRException("Too many steps in routine qtrap");
    }

    public static <T extends Func_Doub_To_Doub> double qsimp(final T func, final double a, final double b) throws NRException {
        return qsimp(func, a, b, 1.0e-10);
    }

    public static <T extends Func_Doub_To_Doub> double qsimp(final T func, final double a, final double b, final double eps)
            throws NRException {
        // Returns the integral of the function or functor func from a to b.
        // The constants EPS can be set to the desired fractional accuracy and
        // JMAX so that 2 to the power JMAX-1 is the maximum allowed number of
        // steps. Integration is performed by Simpson’s rule.
        final int JMAX = 20;
        double s, st, ost = 0.0, os = 0.0;
        Trapzd<T> t = new Trapzd<T>(func, a, b);
        for (int j = 0; j < JMAX; j++) {
            st = t.next();
            s = (4.0 * st - ost) / 3.0; // Compare equation (4.2.4), above.
            if (j > 5) // Avoid spurious early convergence.
                if (abs(s - os) < eps * abs(os) || (s == 0.0 && os == 0.0))
                    return s;
            os = s;
            ost = st;
        }
        throw new NRException("Too many steps in routine qsimp");
    }

    public static <T extends Func_Doub_To_Doub> double qromb(final T func, final double a, final double b) throws NRException {
        return qromb(func, a, b, 1.0e-10);
    }

    public static <T extends Func_Doub_To_Doub> double qromb(final T func, final double a, final double b, final double eps) throws NRException {
        // Returns the integral of the function or functor func from a to b.
        // Integration is performed by Romberg’s method of order 2K, where,
        // e.g., K=2 is Simpson’s rule.
        final int JMAX = 20, JMAXP = JMAX + 1, K = 5;
        // Here EPS is the fractional accuracy desired, as determined by the
        // extrapolation error estimate; JMAX limits the total number of steps;
        // K is the number of points used in the extrapolation.
        final double[] s = doub_vec(JMAX), h = doub_vec(JMAXP); // These store the
                                                              // successive
                                                              // trapezoidal
                                                              // approxi
        Poly_interp polint = new Poly_interp(h, s, K); // mations and their
                                                       // relative stepsizes.
        h[0] = 1.0;
        Trapzd<T> t = new Trapzd<T>(func, a, b);
        for (int j = 1; j <= JMAX; j++) {
            s[j - 1] = t.next();
            if (j >= K) {
                double ss = polint.rawinterp(j - K, 0.0);
                if (abs(polint.dy()) <= eps * abs(ss))
                    return ss;
            }
            h[j] = 0.25 * h[j - 1];
            // This is a key step: The factor is 0.25 even though the stepsize
            // is decreased by only 0.5. This makes the extrapolation a
            // polynomial
            // in h2 as allowed by equation (4.2.1), not just a polynomial in h.
        }
        throw new NRException("Too many steps in routine qromb");
    }

    public static <T extends Func_Doub_To_Doub> double qromo(final Midpnt<T> q) throws NRException {
        return qromo(q, 3.0e-9);
    }

    public static <T extends Func_Doub_To_Doub> double qromo(final Midpnt<T> q, final double eps) throws NRException {
        // Romberg integration on an open interval. Returns the integral of a
        // function using any specified elementary quadrature algorithm q and
        // Romberg’s method. Normally q will be an open formula, not evaluating
        // the function at the endpoints. It is assumed that q triples the
        // number of steps on each call, and that its error series contains
        // only even powers of the number of steps. The routines midpnt,
        // midinf, midsql, midsqu, midexp are possible choices for q. The
        // constants below have the same meanings as in qromb.
        final int JMAX = 14, JMAXP = JMAX + 1, K = 5;
        final double[] h = doub_vec(JMAXP), s = doub_vec(JMAX);
        Poly_interp polint = new Poly_interp(h, s, K);
        h[0] = 1.0;
        for (int j = 1; j <= JMAX; j++) {
            s[j - 1] = q.next();
            if (j >= K) {
                double ss = polint.rawinterp(j - K, 0.0);
                if (abs(polint.dy()) <= eps * abs(ss))
                    return ss;
            }
            h[j] = h[j - 1] / 9.0; // This is where the assumption of step
                                   // tripling and an even
        } // error series is used.
        throw new NRException("Too many steps in routine qromo");
    }

    // Here are the abscissas and weights:
    private static final double x[] = { 0.1488743389816312, 0.4333953941292472, 0.6794095682990244, 0.8650633666889845,
            0.9739065285171717 };
    private static final double w[] = { 0.2955242247147529, 0.2692667193099963, 0.2190863625159821, 0.1494513491505806,
            0.0666713443086881 };

    public static <T extends Func_Doub_To_Doub> double qgaus(final T func, final double a, final double b) throws NRException {
        // Returns the integral of the function or functor func between a and
        // b, by ten-point Gauss- Legendre integration: the function is
        // evaluated exactly ten times at interior points in the range of
        // integration.

        double xm = 0.5 * (b + a);
        double xr = 0.5 * (b - a);
        double s = 0; // Will be twice the average value of the function, since
                      // the
        // ten weights (five numbers above each used twice)
        // sum to 2.
        for (int j = 0; j < 5; j++) {
            double dx = xr * x[j];
            s += w[j] * (func.eval(xm + dx) + func.eval(xm - dx));
        }
        return s *= xr; // Scale the answer to the range of integration.
    }

    public static void gauleg(final double x1, final double x2, final double[] x, final double[] w) {
        // Given the lower and upper limits of integration x1 and x2, this
        // routine returns arrays x[0..n-1] and w[0..n-1] of length n,
        // containing the abscissas and weights of the Gauss-Legendre n-point
        // quadrature formula.

        final double EPS = 1.0e-14; // EPS is the relative precision.
        double z1, z, xm, xl, pp, p3, p2, p1;
        int n = x.length;
        int m = (n + 1) / 2; // The roots are symmetric in the interval, so
        xm = 0.5 * (x2 + x1); // we only have to find half of them.
        xl = 0.5 * (x2 - x1);
        for (int i = 0; i < m; i++) { // Loop over the desired roots.
            z = cos(3.141592654 * (i + 0.75) / (n + 0.5));
            // Starting with this approximation to the ith root, we enter the
            // main loop of refinement
            // by Newton’s method.
            do {
                p1 = 1.0;
                p2 = 0.0;
                for (int j = 0; j < n; j++) { // Loop up the recurrence relation
                                              // to get the
                    p3 = p2; // Legendre polynomial evaluated at z.
                    p2 = p1;
                    p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1);
                }
                // p1 is now the desired Legendre polynomial. We next compute
                // pp, its derivative,
                // by a standard relation involving also p2, the polynomial of
                // one lower order.
                pp = n * (z * p1 - p2) / (z * z - 1.0);
                z1 = z;
                z = z1 - p1 / pp; // Newton’s method.
            } while (abs(z - z1) > EPS);
            x[i] = xm - xl * z; // Scale the root to the desired interval,
            x[n - 1 - i] = xm + xl * z; // and put in its symmetric counterpart.
            w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp); // Compute the weight
            w[n - 1 - i] = w[i]; // and its symmetric counterpart.
        }
    }

    public static void gaulag(final double[] x, final double[] w, final double alf) throws NRException {
        // Given alf, the parameter ? of the Laguerre polynomials, this
        // routine returns arrays x[0..n-1] and w[0..n-1] containing the
        // abscissas and weights of the n-point Gauss-Laguerre quadrature
        // formula. The smallest abscissa is returned in x[0], the largest
        // in x[n-1].

        final int MAXIT = 10;
        final double EPS = 1.0e-14; // EPS is the relative precision.
        int i, its, j;
        double ai, p1, p2 = 0.0, p3, pp = 0.0, z = 0.0, z1;
        int n = x.length;
        for (i = 0; i < n; i++) { // Loop over the desired roots.
            if (i == 0) { // Initial guess for the smallest root.
                z = (1.0 + alf) * (3.0 + 0.92 * alf) / (1.0 + 2.4 * n + 1.8 * alf);
            } else if (i == 1) { // Initial guess for the second root.
                z += (15.0 + 6.25 * alf) / (1.0 + 0.9 * alf + 2.5 * n);
            } else { // Initial guess for the other roots.
                ai = i - 1;
                z += ((1.0 + 2.55 * ai) / (1.9 * ai) + 1.26 * ai * alf / (1.0 + 3.5 * ai)) * (z - x[i - 2])
                        / (1.0 + 0.3 * alf);
            }
            for (its = 0; its < MAXIT; its++) { // Refinement by Newton’s
                                                // method.
                p1 = 1.0;
                p2 = 0.0;
                for (j = 0; j < n; j++) { // Loop up the recurrence relation to
                                          // get the
                    p3 = p2; // Laguerre polynomial evaluated at z.
                    p2 = p1;
                    p1 = ((2 * j + 1 + alf - z) * p2 - (j + alf) * p3) / (j + 1);
                }
                // p1 is now the desired Laguerre polynomial. We next compute
                // pp, its derivative, by a standard relation involving also
                // p2, the polynomial of one lower order.
                pp = (n * p1 - (n + alf) * p2) / z;
                z1 = z;
                z = z1 - p1 / pp; // Newton’s formula.
                if (abs(z - z1) <= EPS)
                    break;
            }
            if (its >= MAXIT)
                throw new NRException("too many iterations in gaulag");
            x[i] = z; // Store the root and the weight.
            w[i] = -exp(gammln(alf + n) - gammln((double) (n))) / (pp * n * p2);
        }
    }

    public static void gauher(final double[] x, final double[] w) throws NRException {
        // This routine returns arrays x[0..n-1] and w[0..n-1] containing the
        // abscissas and weights of // the n-point Gauss-Hermite quadrature
        // formula. The largest abscissa is returned in x[0], the most
        // negative in x[n-1].

        final double EPS = 1.0e-14, PIM4 = 0.7511255444649425;
        // Relative precision and 1= 1=4.
        final int MAXIT = 10; // Maximum iterations.
        int i, its, j, m;
        double p1, p2, p3, pp = 0.0, z = 0.0, z1;
        int n = x.length;
        m = (n + 1) / 2;
        // The roots are symmetric about the origin, so we have to find only
        // half of them.
        for (i = 0; i < m; i++) { // Loop over the desired roots.
            if (i == 0) { // Initial guess for the largest root.
                z = sqrt((double) (2 * n + 1)) - 1.85575 * pow((double) (2 * n + 1), -0.16667);
            } else if (i == 1) { // Initial guess for the second largest root.
                z -= 1.14 * pow((double) (n), 0.426) / z;
            } else if (i == 2) { // Initial guess for the third largest root.
                z = 1.86 * z - 0.86 * x[0];
            } else if (i == 3) { // Initial guess for the fourth largest root.
                z = 1.91 * z - 0.91 * x[1];
            } else { // Initial guess for the other roots.
                z = 2.0 * z - x[i - 2];
            }
            for (its = 0; its < MAXIT; its++) { // Refinement by Newton’s
                                                // method.
                p1 = PIM4;
                p2 = 0.0;
                for (j = 0; j < n; j++) { // Loop up the recurrence relation to
                                          // get
                    // the Hermite polynomial evaluated at
                    // z.
                    p3 = p2;
                    p2 = p1;
                    p1 = z * sqrt(2.0 / (j + 1)) * p2 - sqrt((double) (j) / (j + 1)) * p3;
                }
                // p1 is now the desired Hermite polynomial. We next compute pp,
                // its
                // derivative, by the relation (4.6.21) using p2, the polynomial
                // of one
                // lower order.
                pp = sqrt((double) (2 * n)) * p2;
                z1 = z;
                z = z1 - p1 / pp; // Newton’s formula.
                if (abs(z - z1) <= EPS)
                    break;
            }
            if (its >= MAXIT)
                throw new NRException("too many iterations in gauher");
            x[i] = z; // Store the root
            x[n - 1 - i] = -z; // and its symmetric counterpart.
            w[i] = 2.0 / (pp * pp); // Compute the weight
            w[n - 1 - i] = w[i]; // and its symmetric counterpart.
        }
    }

    public static void gaujac(final double[] x, final double[] w, final double alf, final double bet) throws NRException {
        // Given alf and bet, the parameters ? and ? of the Jacobi polynomials,
        // this routine returns arrays x[0..n-1] and w[0..n-1] containing the
        // abscissas and weights of the n-point Gauss- Jacobi quadrature
        // formula. The largest abscissa is returned in x[0], the smallest
        // in x[n-1].

        final int MAXIT = 10;
        final double EPS = 1.0e-14; // EPS is the relative precision.
        int i, its, j;
        double alfbet, an, bn, r1, r2, r3;
        double a, b, c, p1, p2 = 0.0, p3, pp = 0.0, temp = 0.0, z = 0.0, z1;
        int n = x.length;
        for (i = 0; i < n; i++) { // Loop over the desired roots.
            if (i == 0) { // Initial guess for the largest root.
                an = alf / n;
                bn = bet / n;
                r1 = (1.0 + alf) * (2.78 / (4.0 + n * n) + 0.768 * an / n);
                r2 = 1.0 + 1.48 * an + 0.96 * bn + 0.452 * an * an + 0.83 * an * bn;
                z = 1.0 - r1 / r2;
            } else if (i == 1) { // Initial guess for the second largest root.
                r1 = (4.1 + alf) / ((1.0 + alf) * (1.0 + 0.156 * alf));
                r2 = 1.0 + 0.06 * (n - 8.0) * (1.0 + 0.12 * alf) / n;
                r3 = 1.0 + 0.012 * bet * (1.0 + 0.25 * abs(alf)) / n;
                z -= (1.0 - z) * r1 * r2 * r3;
            } else if (i == 2) { // Initial guess for the third largest root.
                r1 = (1.67 + 0.28 * alf) / (1.0 + 0.37 * alf);
                r2 = 1.0 + 0.22 * (n - 8.0) / n;
                r3 = 1.0 + 8.0 * bet / ((6.28 + bet) * n * n);
                z -= (x[0] - z) * r1 * r2 * r3;
            } else if (i == n - 2) { // Initial guess for the second smallest
                                     // root.
                r1 = (1.0 + 0.235 * bet) / (0.766 + 0.119 * bet);
                r2 = 1.0 / (1.0 + 0.639 * (n - 4.0) / (1.0 + 0.71 * (n - 4.0)));
                r3 = 1.0 / (1.0 + 20.0 * alf / ((7.5 + alf) * n * n));
                z += (z - x[n - 4]) * r1 * r2 * r3;
            } else if (i == n - 1) { // Initial guess for the smallest root.
                r1 = (1.0 + 0.37 * bet) / (1.67 + 0.28 * bet);
                r2 = 1.0 / (1.0 + 0.22 * (n - 8.0) / n);
                r3 = 1.0 / (1.0 + 8.0 * alf / ((6.28 + alf) * n * n));
                z += (z - x[n - 3]) * r1 * r2 * r3;
            } else { // Initial guess for the other roots.
                z = 3.0 * x[i - 1] - 3.0 * x[i - 2] + x[i - 3];
            }
            alfbet = alf + bet;
            for (its = 1; its <= MAXIT; its++) { // Refinement by Newton’s
                                                 // method.
                temp = 2.0 + alfbet; // Start the recurrence with P0 and P1 to
                                     // avoid
                // a division by zero when ? C ? D 0 or 1.
                p1 = (alf - bet + temp * z) / 2.0;
                p2 = 1.0;
                for (j = 2; j <= n; j++) { // Loop up the recurrence relation to
                                           // get the
                    p3 = p2; // Jacobi polynomial evaluated at z.
                    p2 = p1;
                    temp = 2 * j + alfbet;
                    a = 2 * j * (j + alfbet) * (temp - 2.0);
                    b = (temp - 1.0) * (alf * alf - bet * bet + temp * (temp - 2.0) * z);
                    c = 2.0 * (j - 1 + alf) * (j - 1 + bet) * temp;
                    p1 = (b * p2 - c * p3) / a;
                }
                pp = (n * (alf - bet - temp * z) * p1 + 2.0 * (n + alf) * (n + bet) * p2) / (temp * (1.0 - z * z));
                // p1 is now the desired Jacobi polynomial. We next compute pp,
                // its
                // derivative, by a standard relation involving also p2, the
                // polynomial
                // of one lower order.
                z1 = z;
                z = z1 - p1 / pp; // Newton’s formula.
                if (abs(z - z1) <= EPS)
                    break;
            }
            if (its > MAXIT)
                throw new NRException("too many iterations in gaujac");
            x[i] = z; // Store the root and the weight.
            w[i] = exp(gammln(alf + n) + gammln(bet + n) - gammln(n + 1.0) - gammln(n + alfbet + 1.0)) * temp
                    * pow(2.0, alfbet) / (pp * p2);
        }
    }

    public static void gaucof(final double[] a, final double[] b, final double amu0, final double[] x, final double[] w) throws NRException {
        // Computes the abscissas and weights for a Gaussian quadrature formula
        // from the Jacobi matrix. On input, a[0..n-1] and b[0..n-1] are the
        // coefficients of the recurrence relation for the set of monic
        // orthogonal polynomials. The quantity 0 
        // R b
        // a W.x/dx is input as amu0. The abscissas x[0..n-1] are returned
        // in descending order, with the corresponding weights in w[0..n-1].
        // The arrays a and b are modified. Execution can be speeded up by
        // modifying tqli and eigsrt to compute only the zeroth component of
        // each eigenvector.

        int n = a.length;
        for (int i = 0; i < n; i++)
            if (i != 0)
                b[i] = sqrt(b[i]); // Set up superdiagonal of Jacobi
                                        // matrix.
        Symmeig sym = new Symmeig(a, b);
        for (int i = 0; i < n; i++) {
            x[i] = sym.d()[i];
            w[i] = amu0 * sym.z()[0][i] * sym.z()[0][i]; // Equation (4.6.28).
        }
    }

    public static void radau(final double[] a, final double[] b, final double amu0, final double x1, final double[] x, final double[] w)
            throws NRException {
        // Computes the abscissas and weights for a Gauss-Radau quadrature
        // formula. On input, a[0..n-1] and b[0..n-1] are the coefficients of
        // the recurrence relation for the set of monic orthogonal polynomials
        // corresponding to the weight function. (b[0] is not referenced.) The
        // quantity
        // 0 
        // R b
        // a W.x/dx is input as amu0. x1 is input as either endpoint of the
        // interval. The abscissas x[0..n-1] are returned in descending order,
        // with the corresponding weights in w[0..n-1]. The arrays a and b are
        // modified.
        int n = a.length;
        if (n == 1) {
            x[0] = x1;
            w[0] = amu0;
        } else { // Compute pN1 and pN2 by recurrence.
            double p = x1 - a[0];
            double pm1 = 1.0;
            double p1 = p;
            for (int i = 1; i < n - 1; i++) {
                p = (x1 - a[i]) * p1 - b[i] * pm1;
                pm1 = p1;
                p1 = p;
            }
            a[n - 1] = x1 - b[n - 1] * pm1 / p; // Equation (4.6.34).
            gaucof(a, b, amu0, x, w);
        }
    }

    public static void lobatto(final double[] a, final double[] b, final double amu0, final double x1, final double xn, final double[] x,
            final double[] w) throws NRException {
        // Computes the abscissas and weights for a Gauss-Lobatto quadrature
        // formula. On input, the vectors a[0..n-1] and b[0..n-1] are the
        // coefficients of the recurrence relation for the set of monic
        // orthogonal polynomials corresponding to the weight function. (b[0]
        // is not referenced.) The quantity 0 
        // R b
        // a W.x/dx is input as amu0. x1 amd xn are input as the endpoints of
        // the interval. The abscissas x[0..n-1] are returned in descending
        // order, with the corresponding weights in w[0..n-1]. The arrays a
        // and b are modified.
        double det, pl, pr, p1l, p1r, pm1l, pm1r;
        int n = a.length;
        if (n <= 1)
            throw new NRException("n must be bigger than 1 in lobatto");
        pl = x1 - a[0]; // Compute pN1 and pN2 at x1 and xN by recur
        pr = xn - a[0]; // rence.
        pm1l = 1.0;
        pm1r = 1.0;
        p1l = pl;
        p1r = pr;
        for (int i = 1; i < n - 1; i++) {
            pl = (x1 - a[i]) * p1l - b[i] * pm1l;
            pr = (xn - a[i]) * p1r - b[i] * pm1r;
            pm1l = p1l;
            pm1r = p1r;
            p1l = pl;
            p1r = pr;
        }
        det = pl * pm1r - pr * pm1l; // Solve equation (4.6.35).
        a[n - 1] = (x1 * pl * pm1r - xn * pr * pm1l) / det;
        b[n - 1] = (xn - x1) * pl * pr / det;
        gaucof(a, b, amu0, x, w);
    }

    public static <T extends Func_Doub_Doub_Doub_To_Doub> 
            double quad3d(final T func, final double x1, final double x2, 
                            final Func_Doub_To_Doub y1, 
                            final Func_Doub_To_Doub y2, 
                            final Func_Doub_Doub_To_Doub z1, 
                            final Func_Doub_Doub_To_Doub z2) throws NRException {
        // Returns the integral of a user-supplied function func over a
        // three-dimensional region specified by the limits x1, x2, and by
        // the user-supplied functions y1, y2, z1, and z2, as defined
        // in (4.8.2). Integration is performed by calling qgaus recursively.
        NRf1 f1 = new NRf1(y1, y2, z1, z2);
        f1.f2().f3().set_func3d(func);
        return qgaus(f1, x1, x2);
    }
    
}
