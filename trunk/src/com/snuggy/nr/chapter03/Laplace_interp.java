
package com.snuggy.nr.chapter03;

import static com.snuggy.nr.util.Static.*;

import com.snuggy.nr.chapter02.*;
import com.snuggy.nr.refs.*;
import com.snuggy.nr.util.*;

@Deprecated @Broken
public class Laplace_interp extends Linbcg {

    // Object for interpolating missing data in a matrix by solving Laplace’s
    // equation. Call constructor once, then solve one or more times (see text).

    private double[][] mat;
    private int ii, jj;
    private int nn;
    private int[] iter;
    private double[] b, mask;
    private double[] y;

    public Laplace_interp(final double[][] matrix) {
        // Constructor. Values greater than 1.e99 in the input matrix mat are
        // deemed to be missing data. The matrix is not altered until solve
        // is called.
        mat = (matrix);
        ii = (nrows(mat));
        jj = (ncols(mat));
        nn = (ii * jj);
        iter[0] = (0);
        b = doub_arr(nn);
        y = doub_arr(nn);
        mask = doub_arr(nn);
        int i, j, k;
        double vl = 0.;
        for (k = 0; k < nn; k++) { // Fill the r.h.s. vector, the initial guess,
            i = k / jj; // and a mask of the missing data.
            j = k - i * jj;
            if (mat[i][j] < 1.e99) {
                b[k] = y[k] = vl = mat[i][j];
                mask[k] = 1;
            } else {
                b[k] = 0.;
                y[k] = vl;
                mask[k] = 0;
            }
        }
    }

    // void asolve(VecDoub_I &b, VecDoub_O &x, final int  itrnsp);
    // void atimes(VecDoub_I &x, VecDoub_O &r, final int  itrnsp);

    // See definitions below. These are the real algorithmic content.

    public double solve() throws NRException {
        return solve(1.e-6, -1);
    }

    public double solve(final double tol) throws NRException {
        return solve(tol, -1);
    }

    public double solve(final double tol, int itmax) throws NRException {
        // Invoke Linbcg::solve with appropriate arguments. The default
        // argument values will usually work, in which case this routine need
        // be called only once. The original matrix mat is refilled with the
        // interpolated solution.
        int i, j, k;
        double[] err = doub_arr(1);
        if (itmax <= 0)
            itmax = 2 * MAX(ii, jj);
        super.solve(b, y, 1, tol, itmax, iter, err);
        for (k = 0, i = 0; i < ii; i++)
            for (j = 0; j < jj; j++)
                mat[i][j] = y[k++];
        return err[0];
    }

    @Override
    public void asolve(final double[] b, final double[] x, final int itrnsp) {
        // Diagonal preconditioner. (Diagonal elements all unity.)
        int i, n = b.length;
        for (i = 0; i < n; i++)
            x[i] = b[i];
    }

    @Override
    public void atimes(final double[] x, final $double1d r, final int itrnsp) {
        // Sparse matrix, and matrix transpose, multiply. This routine
        // embodies eqs. (3.8.4), (3.8.5), and (3.8.6).
        int i, j, k, n = r.$().length, jjt, it;
        double del;
        for (k = 0; k < n; k++)
            r.$()[k] = 0.;
        for (k = 0; k < n; k++) {
            i = k / jj;
            j = k - i * jj;
            if (mask[k] != 0) { // Measured point, eq. (3.8.5).
                r.$()[k] += x[k];
            } else if (i > 0 && i < ii - 1 && j > 0 && j < jj - 1) { // Interior
                                                                     // point,
                                                                     // eq.
                                                                     // (3.8.4).
                if (itrnsp != 0) {
                    r.$()[k] += x[k];
                    del = -0.25 * x[k];
                    r.$()[k - 1] += del;
                    r.$()[k + 1] += del;
                    r.$()[k - jj] += del;
                    r.$()[k + jj] += del;
                } else {
                    r.$()[k] = x[k] - 0.25 * (x[k - 1] + x[k + 1] + x[k + jj] + x[k - jj]);
                }
            } else if (i > 0 && i < ii - 1) { // Left or right edge, eq.
                                              // (3.8.6).
                if (itrnsp != 0) {
                    r.$()[k] += x[k];
                    del = -0.5 * x[k];
                    r.$()[k - jj] += del;
                    r.$()[k + jj] += del;
                } else {
                    r.$()[k] = x[k] - 0.5 * (x[k + jj] + x[k - jj]);
                }
            } else if (j > 0 && j < jj - 1) { // Top or bottom edge, eq.
                                              // (3.8.6).
                if (itrnsp != 0) {
                    r.$()[k] += x[k];
                    del = -0.5 * x[k];
                    r.$()[k - 1] += del;
                    r.$()[k + 1] += del;
                } else {
                    r.$()[k] = x[k] - 0.5 * (x[k + 1] + x[k - 1]);
                }
            } else { // Corners, eq. (3.8.6).
                jjt = i == 0 ? jj : -jj;
                it = j == 0 ? 1 : -1;
                if (itrnsp != 0) {
                    r.$()[k] += x[k];
                    del = -0.5 * x[k];
                    r.$()[k + jjt] += del;
                    r.$()[k + it] += del;
                } else {
                    r.$()[k] = x[k] - 0.5 * (x[k + jjt] + x[k + it]);
                }
            }
        }
    }
}
