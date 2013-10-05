
package com.snuggy.nr.chapter05;

import static com.snuggy.nr.refs.Refs.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.chapter02.*;
import com.snuggy.nr.refs.*;
import com.snuggy.nr.util.*;

public class Static {

    public static void ddpoly(final double[] c, final double x, final double[] pd) {
        // Given the coefficients of a polynomial of degree nc as an array
        // c[0..nc] of size nc+1 (with c[0] being the constant term), and
        // given a value x, this routine fills an output array pd of size nd+1
        // with the value of the polynomial evaluated at x in pd[0], and the
        // first nd derivatives at x in pd[1..nd].
        int nnd, j, i, nc = c.length - 1, nd = pd.length - 1;
        double cnst = 1.0;
        pd[0] = c[nc];
        for (j = 1; j < nd + 1; j++)
            pd[j] = 0.0;
        for (i = nc - 1; i >= 0; i--) {
            nnd = (nd < (nc - i) ? nd : nc - i);
            for (j = nnd; j > 0; j--)
                pd[j] = pd[j] * x + pd[j - 1];
            pd[0] = pd[0] * x + c[i];
        }
        for (i = 2; i < nd + 1; i++) { // After the first derivative, factorial
                                       // constants come in.
            cnst *= i;
            pd[i] *= cnst;
        }
    }

    public static void poldiv(final double[] u, final double[] v, final $$double1d q, final $$double1d r)
            throws NRException {
        // Divide a polynomial u by a polynomial v, and return the quotient
        // and remainder polynomials in q and r, respectively. The four
        // polynomials are represented as vectors of coefficients, each
        // starting with the constant term. There is no restriction on the
        // relative lengths of u and v, and either may have trailing zeros
        // (represent a lower degree polynomial than its length allows). q and
        // r are returned with the size of u, but will usually have trailing
        // zeros.
        int k, j, n = u.length - 1, nv = v.length - 1;
        while (nv >= 0 && v[nv] == 0.)
            nv--;
        if (nv < 0)
            throw new NRException("poldiv divide by zero polynomial");
        // r = u; May do a resize.
        $$(r, u);
        // q.assign(u.size(), 0.); May do a resize.
        $$(q, new double[u.length]);
        for (k = n - nv; k >= 0; k--) {
            q.$()[k] = r.$()[nv + k] / v[nv];
            for (j = nv + k - 1; j >= k; j--)
                r.$()[j] -= q.$()[k] * v[j - k];
        }
        for (j = nv; j <= n; j++)
            r.$()[j] = 0.0;
    }

    public static <T extends Func_Doub_To_Doub> double dfridr(T func, final double x, final double h,
            final double err_ref[]) throws NRException {
        // Returns the derivative of a function func at a point x by Ridders’
        // method of polynomial extrapolation. The value h is input as an
        // estimated initial stepsize; it need not be small, but rather should
        // be an increment in x over which func changes substantially. An
        // estimate of the error in the derivative is returned as err.
        final int ntab = 10; // Sets maximum size of tableau.
        final double con = 1.4, con2 = (con * con); // Stepsize decreased by CON
                                                    // at each iteration.
        final double big = Double.MAX_VALUE; // numeric_limits<Doub>::max();
        final double safe = 2.0; // Return when error is SAFE worse than the
        int i, j; // best so far.
        double errt, fac, hh, ans = 0.0;
        final double[][] a = doub_mat(ntab, ntab);
        if (h == 0.0)
            throw new NRException("h must be nonzero in dfridr.");
        hh = h;
        a[0][0] = (func.eval(x + hh) - func.eval(x - hh)) / (2.0 * hh);
        err_ref[0] = big;
        for (i = 1; i < ntab; i++) {
            // Successive columns in the Neville tableau will go to smaller
            // stepsizes and higher orders of extrapolation.
            hh /= con;
            // Try new, smaller stepsize.
            a[0][i] = (func.eval(x + hh) - func.eval(x - hh)) / (2.0 * hh);
            fac = con2;
            for (j = 1; j <= i; j++) { // Compute extrapolations of various
                                       // orders, requiring
                // no new function evaluations.
                a[j][i] = (a[j - 1][i] * fac - a[j - 1][i - 1]) / (fac - 1.0);
                fac = con2 * fac;
                errt = MAX(abs(a[j][i] - a[j - 1][i]), abs(a[j][i] - a[j - 1][i - 1]));
                // The error strategy is to compare each new extrapolation to
                // one
                // order lower, both at the present stepsize and the previous
                // one.
                if (errt <= err_ref[0]) { // If error is decreased, save the
                                          // improved answer.
                    err_ref[0] = errt;
                    ans = a[j][i];
                }
            }
            if (abs(a[i][i] - a[i - 1][i - 1]) >= safe * err_ref[0])
                break;
            // If higher order is worse by a significant factor SAFE, then quit
            // early.
        }
        return ans;
    }

    public static void pcshft(final double a, final double b, final double[] d) {
        // Polynomial coefficient shift. Given a coefficient array d[0..n-1],
        // this routine generates a coefficient array g[0..n-1] such that
        // Pn-1
        // kD0 dkyk D
        // Pn-1
        // kD0 gkxk, where x and y are related
        // by (5.8.10), i.e., the interval 1 < y < 1 is mapped to the interval
        // a < x < b. The array g is
        // returned in d.
        int k, j, n = d.length;
        double cnst = 2.0 / (b - a), fac = cnst;
        for (j = 1; j < n; j++) { // First we rescale by the factor const...
            d[j] *= fac;
            fac *= cnst;
        }
        cnst = 0.5 * (a + b); // ...which is then redefined as the desired
                              // shift.
        for (j = 0; j <= n - 2; j++)
            // We accomplish the shift by synthetic division, a miracle
            for (k = n - 2; k >= j; k--)
                // of high-school algebra.
                d[k] -= cnst * d[k + 1];
    }

    public static void ipcshft(final double a, final double b, final double[] d) {
        pcshft(-(2. + b + a) / (b - a), (2. - b - a) / (b - a), d);
    }

    public static Ratfn pade(final double[] cof) throws NRException {
        // Given cof[0..2*n], the leading terms in the power series expansion
        // of a function, solve the linear Pad´e equations to return a Ratfn
        // object that embodies a diagonal rational function approximation to
        // the same function.
        @SuppressWarnings("unused")
        final double BIG = 1.0e99;
        int j, k, n = (cof.length - 1) / 2;
        double sum;
        final double[][] q = doub_mat(n, n);
        @SuppressWarnings("unused")
        final double[][] qlu = doub_mat(n, n);
        @SuppressWarnings("unused")
        final int[] indx = int_arr(n);
        final double[] x = doub_arr(n), y = doub_arr(n), num = doub_arr(n + 1), denom = doub_arr(n + 1);
        for (j = 0; j < n; j++) { // Set up matrix for solving.
            y[j] = cof[n + j + 1];
            for (k = 0; k < n; k++)
                q[j][k] = cof[j - k + n];
        }
        LUdcmp lu = new LUdcmp(q); // Solve by LU decomposition and backsubstitu
        lu.solve(y, x); // tion, with iterative improvement.
        for (j = 0; j < 4; j++)
            lu.mprove(y, x);
        for (k = 0; k < n; k++) { // Calculate the remaining coefficients.
            for (sum = cof[k + 1], j = 0; j <= k; j++)
                sum -= x[j] * cof[k - j];
            y[k] = sum;
        }
        num[0] = cof[0];
        denom[0] = 1.;
        for (j = 0; j < n; j++) { // Copy answers to output.
            num[j + 1] = y[j];
            denom[j + 1] = -x[j];
        }
        return new Ratfn(num, denom);
    }

    public static $$<Ratfn> ratlsq(final Func_Doub_To_Doub fn, final double a, final double b, final int mm,
            final int kk, final double dev_ref[]) throws NRException {
        // Returns a rational function approximation to the function fn in the
        // interval .a; b/. Input quantities mm and kk specify the order of the
        // numerator and denominator, respectively. The maximum absolute
        // deviation of the approximation (insofar as is known) is returned as
        // dev.
        final int NPFAC = 8, MAXIT = 5;
        final double BIG = 1.0e99, PIO2 = 1.570796326794896619;
        int i, it, j, ncof = mm + kk + 1, npt = NPFAC * ncof;
        // Number of points where function is evaluated, i.e., fineness of the
        // mesh.
        double devmax, e, hth, power, sum;
        final double[] bb = doub_arr(npt), coff = doub_arr(ncof), ee = doub_arr(npt), fs = doub_arr(npt), wt = doub_arr(npt), xs = doub_arr(npt);
        final double[][] u = doub_mat(npt, ncof);
        $$<Ratfn> ratbest = $$(new Ratfn(coff, mm + 1, kk + 1));
        dev_ref[0] = BIG;
        for (i = 0; i < npt; i++) { // Fill arrays with mesh abscissas and
                                    // function val
            if (i < (npt / 2) - 1) { // ues.
                hth = PIO2 * i / (npt - 1.0); // At each end, use formula that
                                              // minimizes round
                xs[i] = a + (b - a) * SQR(sin(hth)); // off sensitivity.
            } else {
                hth = PIO2 * (npt - i) / (npt - 1.0);
                xs[i] = b - (b - a) * SQR(sin(hth));
            }
            fs[i] = fn.eval(xs[i]);
            wt[i] = 1.0; // In later iterations we will adjust these weights to
            ee[i] = 1.0; // combat the largest deviations.
        }
        e = 0.0;
        for (it = 0; it < MAXIT; it++) { // Loop over iterations.
            for (i = 0; i < npt; i++) { // Set up the “design matrix” for the
                                        // least-squares
                power = wt[i]; // fit.
                bb[i] = power * (fs[i] + SIGN(e, ee[i]));
                // Key idea here: Fit to fn.x/Ce where the deviation is
                // positive,
                // to fn.x/e where it is negative. Then e is supposed to become
                // an
                // approximation to the equal-ripple deviation.
                for (j = 0; j < mm + 1; j++) {
                    u[i][j] = power;
                    power *= xs[i];
                }
                power = -bb[i];
                for (j = mm + 1; j < ncof; j++) {
                    power *= xs[i];
                    u[i][j] = power;
                }
            }
            SVD svd = new SVD(u); // Singular value decomposition.
            svd.solve(bb, coff);
            // In especially singular or difficult cases, one might here edit
            // the
            // singular values, replacing small values by zero in w[0..ncof-1].
            devmax = sum = 0.0;
            Ratfn rat = new Ratfn(coff, mm + 1, kk + 1);
            for (j = 0; j < npt; j++) { // Tabulate the deviations and revise
                                        // the weights.
                ee[j] = rat.func(xs[j]) - fs[j];
                wt[j] = abs(ee[j]); // Use weighting to emphasize most deviant
                                    // points.
                sum += wt[j];
                if (wt[j] > devmax)
                    devmax = wt[j];
            }
            e = sum / npt; // Update e to be the mean absolute deviation.
            if (devmax <= dev_ref[0]) { // Save only the best coefficient set
                                        // found.
                ratbest.$$(rat);
                dev_ref[0] = devmax;
            }
            // cout << " ratlsq iteration= " << it;
            // cout << " max error= " << setw(10) << devmax << endl;
        }
        return ratbest;
    }
}
