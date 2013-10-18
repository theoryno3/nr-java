
package com.snuggy.nr.chapter09;

import static com.snuggy.nr.util.Static.*;
import com.snuggy.nr.refs.*;
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

    public static <T extends Func_Doub_To_Doub> double rtbis(final T func, final double x1, final double x2, final double xacc)
            throws NRException {
        // Using bisection, return the root of a function or functor func known
        // to lie between x1 and x2. The root will be refined until its accuracy
        // is xacc.
        final int JMAX = 50; // Maximum allowed number of bisections.
        double dx, xmid, rtb;
        double f = func.eval(x1);
        double fmid = func.eval(x2);
        if (f * fmid >= 0.0)
            throw new NRException("Root must be bracketed for bisection in rtbis");
        // rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2); // Orient the search so that f>0
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
}
