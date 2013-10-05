
package com.snuggy.nr.chapter02;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Cholesky {

    // Object for Cholesky decomposition of a matrix A, and related functions.
    private int n;
    private final double[][] el; // Stores the decomposition.

    public Cholesky(final double[][] a) throws NRException {
        n = (nrows(a));
        el = doub_mat(a);
        // Constructor. Given a positive-definite symmetric matrix
        // a[0..n-1][0..n-1], construct and store its Cholesky 
        // decomposition, A D L  LT .
        int i, j, k;
        //final double[] tmp;
        double sum;
        if (ncols(el) != n)
            throw new NRException("need square matrix");
        for (i = 0; i < n; i++) {
            for (j = i; j < n; j++) {
                for (sum = el[i][j], k = i - 1; k >= 0; k--)
                    sum -= el[i][k] * el[j][k];
                if (i == j) {
                    if (sum <= 0.0) // A, with rounding errors, is not
                                    // positive-definite.
                        throw new NRException("Cholesky failed");
                    el[i][i] = sqrt(sum);
                } else
                    el[j][i] = sum / el[i][i];
            }
        }
        for (i = 0; i < n; i++)
            for (j = 0; j < i; j++)
                el[j][i] = 0.;
    }

    public void solve(final double[] b, final double[] x) throws NRException {
        // Solve the set of n linear equations Ax D b, where a is a
        // positive-definite symmetric matrix whose Cholesky
        // decomposition has been stored. b[0..n-1] is input as the
        // right-hand side vector. The solution vector is returned
        // in x[0..n-1].
        int i, k;
        double sum;
        if (b.length != n || x.length != n)
            throw new NRException("bad lengths in Cholesky");
        for (i = 0; i < n; i++) { // Solve L  y D b, storing y in x.
            for (sum = b[i], k = i - 1; k >= 0; k--)
                sum -= el[i][k] * x[k];
            x[i] = sum / el[i][i];
        }
        for (i = n - 1; i >= 0; i--) { // Solve LT  x D y.
            for (sum = x[i], k = i + 1; k < n; k++)
                sum -= el[k][i] * x[k];
            x[i] = sum / el[i][i];
        }
    }

    public void elmult(final double[] y, final double[] b) throws NRException {
        // Multiply L  y D b, where L is the lower triangular matrix
        // in the stored Cholesky decomposition. y[0..n-1] is input.
        // The result is returned in b[0..n-1].
        int i, j;
        if (b.length != n || y.length != n)
            throw new NRException("bad lengths");
        for (i = 0; i < n; i++) {
            b[i] = 0.;
            for (j = 0; j <= i; j++)
                b[i] += el[i][j] * y[j];
        }
    }

    public void elsolve(final double[] b, final double[] y) throws NRException {
        // Solve L  y D b, where L is the lower triangular matrix in
        // the stored Cholesky decomposition. b[0..n-1] is input as
        // the right-hand side vector. The solution vector is returned
        // in y[0..n-1].
        int i, j;
        double sum;
        if (b.length != n || y.length != n)
            throw new NRException("bad lengths");
        for (i = 0; i < n; i++) {
            for (sum = b[i], j = 0; j < i; j++)
                sum -= el[i][j] * y[j];
            y[i] = sum / el[i][i];
        }
    }

    public void inverse(final double[][] ainv_ref[]) {
        // Set ainv[0..n-1][0..n-1] to the matrix inverse of A,
        // the matrix whose Cholesky decomposition has been stored.
        int i, j, k;
        double sum;
        // ainv.resize(n,n);
        ainv_ref[0] = doub_mat(n, n);
        for (i = 0; i < n; i++)
            for (j = 0; j <= i; j++) {
                sum = (i == j ? 1. : 0.);
                for (k = i - 1; k >= j; k--)
                    sum -= el[i][k] * ainv_ref[0][j][k];
                ainv_ref[0][j][i] = sum / el[i][i];
            }
        for (i = n - 1; i >= 0; i--)
            for (j = 0; j <= i; j++) {
                sum = (i < j ? 0. : ainv_ref[0][j][i]);
                for (k = i + 1; k < n; k++)
                    sum -= el[k][i] * ainv_ref[0][j][k];
                ainv_ref[0][i][j] = ainv_ref[0][j][i] = sum / el[i][i];
            }
    }

    public double logdet() {
        // Return the logarithm of the determinant of A, the matrix
        // whose Cholesky decomposition has been stored.
        double sum = 0.;
        for (int i = 0; i < n; i++)
            sum += log(el[i][i]);
        return 2. * sum;
    }

}
