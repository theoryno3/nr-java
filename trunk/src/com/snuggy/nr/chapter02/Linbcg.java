
package com.snuggy.nr.chapter02;

import static com.snuggy.nr.refs.Refs.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.refs.*;
import com.snuggy.nr.util.*;

public abstract class Linbcg {

    // Abstract base class for solving sparse linear equations by the
    // preconditioned biconjugate gradient method. To use, declare a derived
    // class in which the methods atimes and asolve are defined for your
    // problem, along with any data that they need. Then call the solve method.

    public abstract void asolve(final double[] b, final double[] x, final int itrnsp);

    public abstract void atimes(final double[] x, final $$double1d r, final int itrnsp) throws NRException;

    // void solve(VecDoub_I &b, VecDoub_IO &x, final int  itol, final double tol,
    // final int  itmax, Int &iter, Doub &err);

    // Doub snrm(VecDoub_I &sx, final int  itol); Utility used by solve.

    public void solve(final double[] b, final double[] x, final int itol, final double tol, 
                        final int itmax, final $int iter, final $double err) 
            throws NRException {
        // Solves Ax D b for x[0..n-1], given b[0..n-1], by the iterative
        // biconjugate gradient method. On input x[0..n-1] should be set to an
        // initial guess of the solution (or all zeros); itol is 1,2,3, or 4,
        // specifying which convergence test is applied (see text); itmax is
        // the maximum number of allowed iterations; and tol is the desired
        // convergence tolerance. On output, x[0..n-1] is reset to the improved
        // solution, iter is the number of iterations actually taken, and err
        // is the estimated error. The matrix A is referenced only through
        // the user-supplied routines atimes, which computes the product of
        // either A or its transpose on a vector, and asolve, which solves
        // zA  x D b or zA T  x D b for some preconditioner matrix zA
        // (possibly the trivial diagonal part of A). This routine can be
        // called repeatedly, with itmax.n, to monitor how err decreases; or
        // it can be called once with a sufficiently large value of itmax
        // so that convergence to tol is achieved.

        double ak, akden, bk, bkden = 1.0, bknum, bnrm, dxnrm, xnrm, 
               zm1nrm, znrm = 0.0;
        final double EPS = 1.0e-14;
        int j, n = b.length;
        final double[] p = doub_arr(n), pp = doub_arr(n), 
                 rr = doub_arr(n);
        $$double1d r = $$(new double[n]); 
        $$double1d z = $$(new double[n]); 
        $$double1d zz = $$(new double[n]);
        $(iter, 0); // Calculate initial residual.
        atimes(x, r, 0); // Input to atimes is x[0..n-1], output is r[0..n-1];
        // the final 0 indicates that the matrix (not its
        // transpose) is to be used.
        for (j = 0; j < n; j++) {
            r.$()[j] = b[j] - r.$()[j];
            rr[j] = r.$()[j];
        }
        // atimes(r,rr,0); Uncomment this line to get the “minimum resid
        if (itol == 1) { // ual” variant of the algorithm.
            bnrm = snrm(b, itol);
            asolve(r.$(), z.$(), 0); // Input to asolve is r[0..n-1], output is
                             // z[0..n-1];
            // the final 0 indicates that the matrixeA
            // (not its transpose) is to be used.
        } else if (itol == 2) {
            asolve(b, z.$(), 0);
            bnrm = snrm(z.$(), itol);
            asolve(r.$(), z.$(), 0);
        } else if (itol == 3 || itol == 4) {
            asolve(b, z.$(), 0);
            bnrm = snrm(z.$(), itol);
            asolve(r.$(), z.$(), 0);
            znrm = snrm(z.$(), itol);
        } else
            throw new NRException("illegal itol in linbcg");
        while (iter.$() < itmax) { // Main loop.
            //++(iter);
            iter.$(iter.$() + 1);
            asolve(rr, zz.$(), 1); // Final 1 indicates use of transpose matrixeA T
                               // .
            for (bknum = 0.0, j = 0; j < n; j++)
                bknum += z.$()[j] * rr[j];
            // Calculate coefficient bk and direction vectors p and pp.
            if (iter.$() == 1) {
                for (j = 0; j < n; j++) {
                    p[j] = z.$()[j];
                    pp[j] = zz.$()[j];
                }
            } else {
                bk = bknum / bkden;
                for (j = 0; j < n; j++) {
                    p[j] = bk * p[j] + z.$()[j];
                    pp[j] = bk * pp[j] + zz.$()[j];
                }
            }
            bkden = bknum; // Calculate coefficient ak, new iterate x, and new
            atimes(p, z, 0); // residuals r and rr.
            for (akden = 0.0, j = 0; j < n; j++)
                akden += z.$()[j] * pp[j];
            ak = bknum / akden;
            atimes(pp, zz, 1);
            for (j = 0; j < n; j++) {
                x[j] += ak * p[j];
                r.$()[j] -= ak * z.$()[j];
                rr[j] -= ak * zz.$()[j];
            }
            asolve(r.$(), z.$(), 0); // SolveeA  z D r and check stopping criterion.
            if (itol == 1)
                err.$(snrm(r.$(), itol) / bnrm);
            else if (itol == 2)
                err.$(snrm(z.$(), itol) / bnrm);
            else if (itol == 3 || itol == 4) {
                zm1nrm = znrm;
                znrm = snrm(z.$(), itol);
                if (abs(zm1nrm - znrm) > EPS * znrm) {
                    dxnrm = abs(ak) * snrm(p, itol);
                    err.$(znrm / abs(zm1nrm - znrm) * dxnrm);
                } else {
                    err.$(znrm / bnrm); // Error may not be accurate, so loop
                                          // again.
                    continue;
                }
                xnrm = snrm(x, itol);
                if (err.$() <= 0.5 * xnrm)
                    err.$(err.$() / xnrm);
                else {
                    err.$(znrm / bnrm); // Error may not be accurate, so loop
                                          // again.
                    continue;
                }
            }
            if (err.$() <= tol)
                break;
        }
    }

    // The solve routine uses this short utility for computing vector norms:

    private double snrm(final double[] sx, final int itol) {
        // Compute one of two norms for a vector sx[0..n-1], as signaled by
        // itol. Used by solve.
        int i, isamax, n = sx.length;
        double ans;
        if (itol <= 3) {
            ans = 0.0;
            for (i = 0; i < n; i++)
                ans += SQR(sx[i]); // Vector magnitude norm.
            return sqrt(ans);
        } else {
            isamax = 0;
            for (i = 0; i < n; i++) { // Largest component norm.
                if (abs(sx[i]) > abs(sx[isamax]))
                    isamax = i;
            }
            return abs(sx[isamax]);
        }
    }

}
