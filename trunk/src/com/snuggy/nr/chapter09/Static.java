
package com.snuggy.nr.chapter09;

import static com.snuggy.nr.util.Static.*;
import static com.snuggy.nr.util.Complex.*;

import com.snuggy.nr.chapter11.*;
import com.snuggy.nr.refs.*;

import static com.snuggy.nr.refs.Refs.*;

import java.io.*;

import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Static {

    public static <T extends Func_Doub_To_Doub> double zbrent(final T func, final double x1, final double x2,
            final double tol) throws NRException {
        // Using Brent’s method, return the root of a function or functor func
        // known to lie between x1 and x2. The root will be refined until its
        // accuracy is tol.
        final int ITMAX = 100; // Maximum allowed number of iterations.
        final double EPS = EPS(); // numeric_limits<Doub>::epsilon();
        // Machine floating-point precision.
        double a = x1, b = x2, c = x2, d = 0.0, e = 0.0, fa = func.eval(a), fb = func.eval(b), fc, p, q, r, s, tol1, xm;
        if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
            throw new NRException("Root must be bracketed in zbrent");
        fc = fb;
        for (int iter = 0; iter < ITMAX; iter++) {
            if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
                c = a; // Rename a, b, c and adjust bounding interval
                fc = fa; // d.
                e = d = b - a;
            }
            if (abs(fc) < abs(fb)) {
                a = b;
                b = c;
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }
            tol1 = 2.0 * EPS * abs(b) + 0.5 * tol; // Convergence check.
            xm = 0.5 * (c - b);
            if (abs(xm) <= tol1 || fb == 0.0)
                return b;
            if (abs(e) >= tol1 && abs(fa) > abs(fb)) {
                s = fb / fa; // Attempt inverse quadratic interpolation.
                if (a == c) {
                    p = 2.0 * xm * s;
                    q = 1.0 - s;
                } else {
                    q = fa / fc;
                    r = fb / fc;
                    p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                    q = (q - 1.0) * (r - 1.0) * (s - 1.0);
                }
                if (p > 0.0)
                    q = -q; // Check whether in bounds.
                p = abs(p);
                double min1 = 3.0 * xm * q - abs(tol1 * q);
                double min2 = abs(e * q);
                if (2.0 * p < (min1 < min2 ? min1 : min2)) {
                    e = d; // Accept interpolation.
                    d = p / q;
                } else {
                    d = xm; // Interpolation failed, use bisection.
                    e = d;
                }
            } else { // Bounds decreasing too slowly, use bisection.
                d = xm;
                e = d;
            }
            a = b; // Move last best guess to a.
            fa = fb;
            if (abs(d) > tol1) // Evaluate new trial root.
                b += d;
            else
                b += SIGN(tol1, xm);
            fb = func.eval(b);
        }
        throw new NRException("Maximum number of iterations exceeded in zbrent");
    }

    public static <T extends Func_Doub_To_Doub> void scrsho(final T fx) throws NRException, NumberFormatException,
            IOException {
        // Graph the function or functor fx over the prompted-for interval
        // x1,x2. Query for another plot until the user signals satisfaction.
        final int RES = 500; // Number of function evaluations for each plot.
        final double XLL = 75., XUR = 525., YLL = 250., YUR = 700.; // Corners
                                                                    // of plot,
                                                                    // in
                                                                    // points.
        String plotfilename = null; // tmpnam(NULL);
        double[] xx = doub_vec(RES), yy = doub_vec(RES);
        double x1, x2;
        int i;
        for (;;) {
            double ymax = -9.99e99, ymin = 9.99e99, del;
            // cout << endl << "Enter x1 x2 (x1=x2 to stop):" << endl;
            System.out.println();
            System.out.println("Enter x1 x2 (x1=x2 to stop):");
            // cin >> x1 >> x2; // Query for another plot, quit if x1=x2.
            BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
            x1 = Double.valueOf(br.readLine());
            x2 = Double.valueOf(br.readLine());
            if (x1 == x2)
                break;
            for (i = 0; i < RES; i++) { // Evaluate the function at equal
                                        // intervals. Find
                xx[i] = x1 + i * (x2 - x1) / (RES - 1.); // the largest and
                                                         // smallest values.
                yy[i] = fx.eval(xx[i]);
                if (yy[i] > ymax)
                    ymax = yy[i];
                if (yy[i] < ymin)
                    ymin = yy[i];
            }
            del = 0.05 * ((ymax - ymin) + (ymax == ymin ? abs(ymax) : 0.));
            // Plot commands, following, are in PSplot syntax (22.1). You can
            // substitute commands for your favorite plotting package.
            PSpage pg = new PSpage(plotfilename);
            PSplot plot = new PSplot(pg, XLL, XUR, YLL, YUR);
            plot.setlimits(x1, x2, ymin - del, ymax + del);
            plot.frame();
            plot.autoscales();
            plot.lineplot(xx, yy);
            if (ymax * ymin < 0.)
                plot.lineseg(x1, 0., x2, 0.);
            plot.display();
        }
        // remove(plotfilename);
        (new File(plotfilename)).delete();
    }

    static class PSpage {
        public PSpage(String s) {
        }
    }

    static class PSplot {
        public PSplot(PSpage pg, double xLL, double xUR, double yLL, double yUR) {
        }

        public void display() {
        }

        public void lineseg(double x1, double d, double x2, double e) {
        }

        public void lineplot(double[] xx, double[] yy) {
        }

        public void autoscales() {
        }

        public void frame() {
        }

        public void setlimits(double x1, double x2, double d, double e) {
        }
    }

    public static <T extends Func_Doub_To_Doub> boolean zbrac(final T func, final $double x1, final $double x2)
            throws NRException {
        // Given a function or functor func and an initial guessed range x1 to
        // x2, the routine expands the range geometrically until a root is
        // bracketed by the returned values x1 and x2 (in which case zbrac
        // returns true) or until the range becomes unacceptably large (in which
        // case zbrac returns false).
        final int NTRY = 50;
        final double FACTOR = 1.6;
        if (x1 == x2)
            throw new NRException("Bad initial range in zbrac");
        double f1 = func.eval(x1.$());
        double f2 = func.eval(x2.$());
        for (int j = 0; j < NTRY; j++) {
            if (f1 * f2 < 0.0)
                return true;
            if (abs(f1) < abs(f2)) {
                x1.$(x1.$() + FACTOR * (x1.$() - x2.$()));
                f1 = func.eval(x1.$());
            } else {
                x2.$(x2.$() + FACTOR * (x2.$() - x1.$()));
                f2 = func.eval(x2.$());
            }
        }
        return false;
    }

    public static <T extends Func_Doub_To_Doub> void zbrak(final T fx, final double x1, final double x2, final int n,
            final $double1d xb1, final $double1d xb2, final $int nroot) throws NRException {
        // Given a function or functor fx defined on the interval [x1,x2],
        // subdivide the interval into n equally spaced segments, and search
        // for zero crossings of the function. nroot will be set to the number
        // of bracketing pairs found. If it is positive, the arrays
        // xb1[0..nroot-1] and xb2[0..nroot-1] will be filled sequentially with
        // any bracketing pairs that are found. On input, these vectors may have
        // any size, including zero; they will be resized to nroot.
        int nb = 20;
        // xb1.resize(nb);
        xb1.$(doub_vec(nb));
        // xb2.resize(nb);
        xb2.$(doub_vec(nb));
        nroot.$(0);
        double dx = (x2 - x1) / n; // Determine the spacing appropriate to the
                                   // mesh.
        double x = x1;
        double fp = fx.eval(x1);
        for (int i = 0; i < n; i++) { // Loop over all intervals
            double fc = fx.eval(x += dx);
            if (fc * fp <= 0.0) { // If a sign change occurs, then record values
                                  // for the
                xb1.$()[nroot.$()] = x - dx; // bounds.
                xb2.$()[nroot.$()] = x;
                nroot.$(nroot.$() + 1);
                if (nroot.$() == nb) {
                    double[] tempvec1 = doub_vec(xb1.$()), tempvec2 = doub_vec(xb2.$());
                    // xb1.resize(2 * nb);
                    xb1.$(doub_vec(2 * nb));
                    // xb2.resize(2 * nb);
                    xb2.$(doub_vec(2 * nb));
                    for (int j = 0; j < nb; j++) {
                        xb1.$()[j] = tempvec1[j];
                        xb2.$()[j] = tempvec2[j];
                    }
                    nb *= 2;
                }
            }
            fp = fc;
        }
    }

    public static <T extends Func_Doub_To_Doub> double rtbis(final T func, final double x1, final double x2,
            final double xacc) throws NRException {
        // Using bisection, return the root of a function or functor func known
        // to lie between x1 and x2. The root will be refined until its accuracy
        // is xacc.
        final int JMAX = 50; // Maximum allowed number of bisections.
        double dx, xmid, rtb;
        double f = func.eval(x1);
        double fmid = func.eval(x2);
        if (f * fmid >= 0.0)
            throw new NRException("Root must be bracketed for bisection in rtbis");
        // rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2); // Orient the search
        // so that f>0
        if (f < 0.0) {
            dx = x2 - x1;
            rtb = (x1);
        } else {
            dx = x1 - x2;
            rtb = (x2);
        }
        for (int j = 0; j < JMAX; j++) { // lies at x+dx.
            fmid = func.eval(xmid = rtb + (dx *= 0.5)); // Bisection loop.
            if (fmid <= 0.0)
                rtb = xmid;
            if (abs(dx) < xacc || fmid == 0.0)
                return rtb;
        }
        throw new NRException("Too many bisections in rtbis");
    }

    public static <T extends Func_Doub_To_Doub> double rtflsp(final T func, final double x1, final double x2,
            final double xacc) throws NRException {
        // Using the false-position method, return the root of a function or
        // functor func known to lie between x1 and x2. The root is refined
        // until its accuracy is xacc.
        final int MAXIT = 30; // Set to the maximum allowed number of
                              // iterations.
        double xl, xh, del;
        double fl = func.eval(x1);
        double fh = func.eval(x2); // Be sure the interval brackets a root.
        if (fl * fh > 0.0)
            throw new NRException("Root must be bracketed in rtflsp");
        if (fl < 0.0) { // Identify the limits so that xl corresponds to the low
            xl = x1; // side.
            xh = x2;
        } else {
            xl = x2;
            xh = x1;
            // SWAP(fl, fh);
            double temp = fl;
            fl = fh;
            fh = temp;
        }
        double dx = xh - xl;
        for (int j = 0; j < MAXIT; j++) { // False-position loop.
            double rtf = xl + dx * fl / (fl - fh); // Increment with respect to
                                                   // latest value.
            double f = func.eval(rtf);
            if (f < 0.0) { // Replace appropriate limit.
                del = xl - rtf;
                xl = rtf;
                fl = f;
            } else {
                del = xh - rtf;
                xh = rtf;
                fh = f;
            }
            dx = xh - xl;
            if (abs(del) < xacc || f == 0.0)
                return rtf; // Convergence.
        }
        throw new NRException("Maximum number of iterations exceeded in rtflsp");
    }

    public static <T extends Func_Doub_To_Doub> double rtsec(final T func, final double x1, final double x2,
            final double xacc) throws NRException {
        // Using the secant method, return the root of a function or functor
        // func thought to lie between x1 and x2. The root is refined until its
        // accuracy is xacc.
        final int MAXIT = 30; // Maximum allowed number of iterations.
        double xl, rts;
        double fl = func.eval(x1);
        double f = func.eval(x2);
        if (abs(fl) < abs(f)) { // Pick the bound with the smaller function
                                // value as
            rts = x1; // the most recent guess.
            xl = x2;
            // SWAP(fl, f);
            double temp = fl;
            fl = f;
            f = temp;
        } else {
            xl = x1;
            rts = x2;
        }
        for (int j = 0; j < MAXIT; j++) { // Secant loop.
            double dx = (xl - rts) * f / (f - fl); // Increment with respect to
                                                   // latest value.
            xl = rts;
            fl = f;
            rts += dx;
            f = func.eval(rts);
            if (abs(dx) < xacc || f == 0.0)
                return rts; // Convergence.
        }
        throw new NRException("Maximum number of iterations exceeded in rtsec");
    }

    public static <T extends Func_Doub_To_Doub> double zriddr(final T func, final double x1, final double x2,
            final double xacc) throws NRException {
        // Using Ridders’ method, return the root of a function or functor func
        // known to lie between x1 and x2. The root will be refined to an
        // approximate accuracy xacc.
        final int MAXIT = 60;
        double fl = func.eval(x1);
        double fh = func.eval(x2);
        if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
            double xl = x1;
            double xh = x2;
            double ans = -9.99e99; // Any highly unlikely value, to simplify
                                   // logic
            for (int j = 0; j < MAXIT; j++) { // below.
                double xm = 0.5 * (xl + xh);
                double fm = func.eval(xm); // First of two function evaluations
                                           // per it
                double s = sqrt(fm * fm - fl * fh); // eration.
                if (s == 0.0)
                    return ans;
                double xnew = xm + (xm - xl) * ((fl >= fh ? 1.0 : -1.0) * fm / s); // Updating
                                                                                   // formula.
                if (abs(xnew - ans) <= xacc)
                    return ans;
                ans = xnew;
                double fnew = func.eval(ans); // Second of two function
                                              // evaluations per
                if (fnew == 0.0)
                    return ans; // iteration.
                if (SIGN(fm, fnew) != fm) { // Bookkeeping to keep the root
                                            // bracketed
                    xl = xm; // on next iteration.
                    fl = fm;
                    xh = ans;
                    fh = fnew;
                } else if (SIGN(fl, fnew) != fl) {
                    xh = ans;
                    fh = fnew;
                } else if (SIGN(fh, fnew) != fh) {
                    xl = ans;
                    fl = fnew;
                } else
                    throw new NRException("never get here.");
                if (abs(xh - xl) <= xacc)
                    return ans;
            }
            throw new NRException("zriddr exceed maximum iterations");
        } else {
            if (fl == 0.0)
                return x1;
            if (fh == 0.0)
                return x2;
            throw new NRException("root must be bracketed in zriddr.");
        }
    }

    /*
     * public static <T extends Func_Doub_To_Doub> double zbrent(final T func,
     * final double x1, final double x2, final double tol) { // Using Brent’s
     * method, return the root of a function or functor func // known to lie
     * between x1 and x2. The root will be refined until its // accuracy is tol.
     * final int ITMAX = 100; // Maximum allowed number of iterations. final
     * double EPS = EPS(); // numeric_limits<double>::epsilon(); // Machine
     * floating-point precision. double a = x1, b = x2, c = x2, d, e, fa =
     * func.eval(a), fb = func.eval(b), fc, p, q, r, s, tol1, xm; if ((fa > 0.0
     * && fb > 0.0) || (fa < 0.0 && fb < 0.0)) throw new
     * NRException("Root must be bracketed in zbrent"); fc = fb; for (int iter =
     * 0; iter < ITMAX; iter++) { if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc
     * < 0.0)) { c = a; // Rename a, b, c and adjust bounding interval fc = fa;
     * // d. e = d = b - a; } if (abs(fc) < abs(fb)) { a = b; b = c; c = a; fa =
     * fb; fb = fc; fc = fa; } tol1 = 2.0 * EPS * abs(b) + 0.5 * tol; //
     * Convergence check. xm = 0.5 * (c - b); if (abs(xm) <= tol1 || fb == 0.0)
     * return b; if (abs(e) >= tol1 && abs(fa) > abs(fb)) { s = fb / fa; //
     * Attempt inverse quadratic interpolation. if (a == c) { p = 2.0 * xm * s;
     * q = 1.0 - s; } else { q = fa / fc; r = fb / fc; p = s * (2.0 * xm * q *
     * (q - r) - (b - a) * (r - 1.0)); q = (q - 1.0) * (r - 1.0) * (s - 1.0); }
     * if (p > 0.0) q = -q; // Check whether in bounds. p = abs(p); double min1
     * = 3.0 * xm * q - abs(tol1 * q); double min2 = abs(e * q); if (2.0 * p <
     * (min1 < min2 ? min1 : min2)) { e = d; // Accept interpolation. d = p / q;
     * } else { d = xm; // Interpolation failed, use bisection. e = d; } } else
     * { // Bounds decreasing too slowly, use bisection. d = xm; e = d; } a = b;
     * // Move last best guess to a. fa = fb; if (abs(d) > tol1) // Evaluate new
     * trial root. b += d; else b += SIGN(tol1, xm); fb = func.eval(b); } throw
     * new NRException("Maximum number of iterations exceeded in zbrent"); }
     */

    public static <T extends Funcd> double rtnewt(final T funcd, final double x1, final double x2, final double xacc)
            throws NRException {
        // Using the Newton-Raphson method, return the root of a function known
        // to lie in the interval OEx1; x2 . The root will be refined until its
        // accuracy is known within xacc. funcd is a usersupplied struct that
        // returns the function value as a functor and the first derivative of
        // the function at the point x as the function df (see text).
        final int JMAX = 20; // Set to maximum number of iterations.
        double rtn = 0.5 * (x1 + x2); // Initial guess.
        for (int j = 0; j < JMAX; j++) {
            double f = funcd.eval(rtn);
            double df = funcd.df(rtn);
            double dx = f / df;
            rtn -= dx;
            if ((x1 - rtn) * (rtn - x2) < 0.0)
                throw new NRException("Jumped out of brackets in rtnewt");
            if (abs(dx) < xacc)
                return rtn; // Convergence.
        }
        throw new NRException("Maximum number of iterations exceeded in rtnewt");
    }

    public static <T extends Funcd> double rtsafe(final T funcd, final double x1, final double x2, final double xacc)
            throws NRException {
        // Using a combination of Newton-Raphson and bisection, return the root
        // of a function bracketed between x1 and x2. The root will be refined
        // until its accuracy is known within xacc. funcd is a user-supplied
        // struct that returns the function value as a functor and the first
        // derivative of the function at the point x as the function df (see
        // text).
        final int MAXIT = 100; // Maximum allowed number of iterations.
        double xh, xl;
        double fl = funcd.eval(x1);
        double fh = funcd.eval(x2);
        if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
            throw new NRException("Root must be bracketed in rtsafe");
        if (fl == 0.0)
            return x1;
        if (fh == 0.0)
            return x2;
        if (fl < 0.0) { // Orient the search so that f.xl/ < 0.
            xl = x1;
            xh = x2;
        } else {
            xh = x1;
            xl = x2;
        }
        double rts = 0.5 * (x1 + x2); // Initialize the guess for root,
        double dxold = abs(x2 - x1); // the “stepsize before last,”
        double dx = dxold; // and the last step.
        double f = funcd.eval(rts);
        double df = funcd.df(rts);
        for (int j = 0; j < MAXIT; j++) { // Loop over allowed iterations.
            if ((((rts - xh) * df - f) * ((rts - xl) * df - f) > 0.0) // Bisect
                                                                      // if
                                                                      // Newton
                                                                      // out of
                                                                      // range,
                    || (abs(2.0 * f) > abs(dxold * df))) { // or not decreasing
                                                           // fast enough.
                dxold = dx;
                dx = 0.5 * (xh - xl);
                rts = xl + dx;
                if (xl == rts)
                    return rts; // Change in root is negligible.
            } else { // Newton step acceptable. Take it.
                dxold = dx;
                dx = f / df;
                double temp = rts;
                rts -= dx;
                if (temp == rts)
                    return rts;
            }
            if (abs(dx) < xacc)
                return rts; // Convergence criterion.
            double f_ = funcd.eval(rts);
            @SuppressWarnings("unused")
            double df_ = funcd.df(rts);
            // The one new function evaluation per iteration.
            if (f_ < 0.0) // Maintain the bracket on the root.
                xl = rts;
            else
                xh = rts;
        }
        throw new NRException("Maximum number of iterations exceeded in rtsafe");
    }

    // static final double frac[MR+1]=
    static final double frac[] = { 0.0, 0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0 };

    public static void laguer(final Complex[] a, final $<Complex> x, final $int its) throws NRException {
        // Given the m+1 complex coefficients a[0..m] of the polynomial
        // Pmi D0 aOEi xi , and given a complex value x, this routine improves
        // x by Laguerre’s method until it converges, within the achievable
        // roundoff limit, to a root of the given polynomial. The number of
        // iterations taken is returned as its.
        final int MR = 8, MT = 10, MAXIT = MT * MR;
        final double EPS = EPS(); // numeric_limits<double>::epsilon();
        // Here EPS is the estimated fractional roundoff error. We try to break
        // (rare) limit cycles with MR different fractional values, once every
        // MT steps, for MAXIT total allowed iterations.
        // Fractions used to break a limit cycle.
        Complex dx, x1, b, d, f, g, h, sq, gp, gm, g2;
        int m = a.length - 1;
        for (int iter = 1; iter <= MAXIT; iter++) { // Loop over iterations up
                                                    // to allowed maximum.
            its.$(iter);
            b = a[m];
            double err = abs(b);
            // d=f=0.0;
            f = complex(0.0);
            d = complex(f);
            double abx = abs(x.$());
            for (int j = m - 1; j >= 0; j--) { // Efficient computation of the
                                               // polynomial and
                // f=x*f+d; // its first two derivatives. f stores P00=2.
                f = plus(times(x.$(), f), d);
                // d=x*d+b;
                d = plus(times(x.$(), d), b);
                // b=x*b+a[j];
                b = plus(times(x.$(), b), a[j]);
                err = abs(b) + abx * err;
            }
            err *= EPS;
            // Estimate of roundoff error in evaluating polynomial.
            if (abs(b) <= err)
                return; // We are on the root.
            // g=d/b; // The generic case: Use Laguerre’s formula.
            g = divide(d, b);
            // g2=g*g;
            g2 = times(g, g);
            // h=g2-2.0*f/b;
            h = minus(g2, times(2.0, divide(f, b)));
            // sq=sqrt(Doub(m-1)*(Doub(m)*h-g2));
            sq = sqrt(times(Doub(m - 1), minus(times(Doub(m), h), g2)));
            // gp=g+sq;
            gp = plus(g, sq);
            // gm=g-sq;
            gm = minus(g, sq);
            double abp = abs(gp);
            double abm = abs(gm);
            if (abp < abm)
                gp = gm;
            dx = MAX(abp, abm) > 0.0 ? divide(Doub(m), gp) : polar(1 + abx, Doub(iter));
            // x1=x-dx;
            x1 = minus(x.$(), dx);
            if (x == x1)
                return; // Converged.
            if (iter % MT != 0)
                x.$(x1);
            else
                x.$(minus(x.$(), times(frac[iter / MT], dx)));
            // Every so often we take a fractional step, to break any limit
            // cycle
            // (itself a rare occurrence).
        }
        throw new NRException("too many iterations in laguer");
        // Very unusual; can occur only for complex roots. Try a different
        // starting guess.
    }

    public static void zroots(final Complex[] a, final Complex[] roots, final $boolean polish) throws NRException,
            InstantiationException, IllegalAccessException {
        // Given the m+1 complex coefficients a[0..m] of the polynomial
        // Pm
        // iD0 a.i/xi , this routine successively calls laguer and finds all
        // m complex roots in roots[0..m-1]. The boolean variable polish
        // should be input as true if polishing (also by Laguerre’s method) is
        // desired, false if the roots will be subsequently polished by other
        // means.
        final double EPS = 1.0e-14; // A small number.
        $int its = $(0);
        int i;
        $<Complex> x = $(complex(0.0));
        Complex b, c;
        int m = a.length - 1;
        Complex[] ad = obj_vec(Complex.class, m + 1);
        for (int j = 0; j <= m; j++)
            ad[j] = a[j]; // Copy of coefficients for successive deflation.
        for (int j = m - 1; j >= 0; j--) { // Loop over each root to be found.
            x = $(complex(0.0)); // Start at zero to favor convergence to small-
            Complex[] ad_v = obj_vec(Complex.class, j + 2); // est remaining
                                                            // root, and return
                                                            // the root.
            for (int jj = 0; jj < j + 2; jj++)
                ad_v[jj] = ad[jj];
            laguer(ad_v, x, its);
            if (abs(imag(x.$())) <= 2.0 * EPS * abs(real(x.$())))
                x.$(complex(real(x.$()), 0.0));
            roots[j] = complex(x.$());
            b = ad[j + 1]; // Forward deflation.
            for (int jj = j; jj >= 0; jj--) {
                c = ad[jj];
                ad[jj] = b;
                // b=x*b+c;
                b = plus(times(x.$(), b), c);
            }
        }
        if (polish.$())
            for (int j = 0; j < m; j++)
                // Polish the roots using the undeflated coeffi
                // laguer(a, roots[j], its); // cients.
                laguer(a, $_(roots, j), its); // cients.
        for (int j = 1; j < m; j++) { // Sort roots by their real parts by
                                      // straight in
            x.$(roots[j]); // sertion.
            for (i = j - 1; i >= 0; i--) {
                if (real(roots[i]) <= real(x.$()))
                    break;
                roots[i + 1] = roots[i];
            }
            roots[i + 1] = complex(x.$());
        }
    }

    public static void zrhqr(final double[] a, final Complex[] rt) throws InstantiationException,
            IllegalAccessException, NRException {
        // Find all the roots of a polynomial with real coefficients,
        // Pm
        // iD0 a.i/xi , given the coefficients
        // a[0..m]. The method is to construct an upper Hessenberg matrix whose
        // eigenvalues are the desired roots and then use the routine Unsymmeig.
        // The roots are returned in the complex vector rt[0..m-1], sorted in
        // descending order by their real parts.
        int m = a.length - 1;
        double[][] hess = doub_mat(m, m);
        for (int k = 0; k < m; k++) { // Construct the matrix.
            hess[0][k] = -a[m - k - 1] / a[m];
            for (int j = 1; j < m; j++)
                hess[j][k] = 0.0;
            if (k != m - 1)
                hess[k + 1][k] = 1.0;
        }
        Unsymmeig h = new Unsymmeig(hess, false, true); // Find its
                                                        // eigenvalues.
        for (int j = 0; j < m; j++)
            rt[j] = complex(h.wri(j));
    }

}
