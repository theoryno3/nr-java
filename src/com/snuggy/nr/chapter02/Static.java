
package com.snuggy.nr.chapter02;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Static {

    public static void gaussj(final double[][] a, final double[][] b) throws NRException {
        // Linear equation solution by Gauss-Jordan elimination,
        // equation (2.1.1) above. The input matrix is a[0..n-1][0..n-1].
        // b[0..n-1][0..m-1] is input containing the m right-hand side vectors.
        // On output, a is replaced by its matrix inverse, and b is replaced by
        // the corresponding set of solution vectors.
        int i, icol = 0, irow = 0, j, k, l, ll, n = nrows(a), m = ncols(b);
        double big, dum, pivinv;
        final int[] indxc = int_arr(n);
        final int[] indxr = int_arr(n);
        final int[] ipiv = int_arr(n); // These integer arrays are used
                                 // for bookkeeping on the pivoting.
        for (j = 0; j < n; j++)
            ipiv[j] = 0;
        for (i = 0; i < n; i++) { // This is the main loop over the columns to
                                  // be
            big = 0.0; // reduced.
            for (j = 0; j < n; j++)
                // This is the outer loop of the search for a pivot
                if (ipiv[j] != 1) // element.
                    for (k = 0; k < n; k++) {
                        if (ipiv[k] == 0) {
                            if (abs(a[j][k]) >= big) {
                                big = abs(a[j][k]);
                                irow = j;
                                icol = k;
                            }
                        }
                    }
            ++(ipiv[icol]);
            // We now have the pivot element, so we interchange rows, if needed,
            // to
            // put the pivot element on the diagonal. The columns are not
            // physically
            // interchanged, only relabeled: indxc[i], the column of the .iC1/th
            // pivot
            // element, is the .iC1/th column that is reduced, while indxr[i] is
            // the
            // row in which that pivot element was originally located.
            // If indxr[i] ¤ indxc[i], there is an implied column interchange.
            // With
            // this form of bookkeeping, the solution b’s will end up in the
            // correct
            // order, and the inverse matrix will be scrambled by columns.
            if (irow != icol) {
                for (l = 0; l < n; l++)
                    SWAP(a, irow, l, icol, l);
                for (l = 0; l < m; l++)
                    SWAP(b, irow, l, icol, l);
            }
            indxr[i] = irow; // We are now ready to divide the pivot row by the
            // pivot element, located at irow and icol.
            indxc[i] = icol;
            if (a[icol][icol] == 0.0)
                throw new NRException("gaussj: Singular Matrix");
            pivinv = 1.0 / a[icol][icol];
            a[icol][icol] = 1.0;
            for (l = 0; l < n; l++)
                a[icol][l] *= pivinv;
            for (l = 0; l < m; l++)
                b[icol][l] *= pivinv;
            for (ll = 0; ll < n; ll++)
                // Next, we reduce the rows...
                if (ll != icol) { // ...except for the pivot one, of course.
                    dum = a[ll][icol];
                    a[ll][icol] = 0.0;
                    for (l = 0; l < n; l++)
                        a[ll][l] -= a[icol][l] * dum;
                    for (l = 0; l < m; l++)
                        b[ll][l] -= b[icol][l] * dum;
                }
        }
        // This is the end of the main loop over columns of the reduction. It
        // only remains to unscramble the solution in view of the column
        // interchanges. We do this by interchanging pairs of columns in the
        // reverse order that the permutation was built up.
        for (l = n - 1; l >= 0; l--) {
            if (indxr[l] != indxc[l])
                for (k = 0; k < n; k++)
                    SWAP(a, k, indxr[l], k, indxc[l]);
        } // And we are done.
    }

    public static void gaussj(final double[][] a) throws NRException {
        // Overloaded version with no right-hand sides. Replaces a by its
        // inverse.
        final double[][] b = doub_mat(nrows(a), 0); // Dummy vector with zero
                                                // columns.
        gaussj(a, b);
    }
    
    public static void tridag(final double[] a, final double[] b, final double[] c, final double[] r, final double[] u) throws NRException {
        // Solves for a vector u[0..n-1] the tridiagonal linear set given by
        // equation (2.4.1). a[0..n-1], b[0..n-1], c[0..n-1], and r[0..n-1] are
        // input vectors and are not modified.
        int j, n = a.length;
        double bet;
        final double[] gam = doub_arr(n); // One vector of workspace, gam, is
                                      // needed.
        if (b[0] == 0.0)
            throw new NRException("Error 1 in tridag");
        // If this happens, then you should rewrite your equations as a
        // set of order N  1, with u1 trivially eliminated.
        u[0] = r[0] / (bet = b[0]);
        for (j = 1; j < n; j++) { // Decomposition and forward substitution.
            gam[j] = c[j - 1] / bet;
            bet = b[j] - a[j] * gam[j];
            if (bet == 0.0)
                throw new NRException("Error 2 in tridag"); // Algorithm fails;
                                                            // see below.
            u[j] = (r[j] - a[j] * u[j - 1]) / bet;
        }
        for (j = (n - 2); j >= 0; j--)
            u[j] -= gam[j + 1] * u[j + 1]; // Backsubstitution.
    }

    public static void banmul(final double[][] a, final int m1, final int m2, final double[] x, final double[] b) {
        // Matrix multiply b D A  x, where A is band-diagonal with m1 rows
        // below the diagonal and m2 rows above. The input vector is x[0..n-1]
        // and the output vector is b[0..n-1]. The array a[0..n-1][0..m1+m2]
        // stores A as follows: The diagonal elements are in a[0..n-1][m1].
        // Subdiagonal elements are in a[j ..n-1][0..m1-1] with j > 0
        // appropriate to the number of elements on each subdiagonal.
        // Superdiagonal elements are in a[0..j ][m1+1..m1+m2] with
        // j < n-1 appropriate to the number of elements on each superdiagonal.
        int i, j, k, tmploop, n = nrows(a);
        for (i = 0; i < n; i++) {
            k = i - m1;
            tmploop = MIN(m1 + m2 + 1, (int) (n - k));
            b[i] = 0.0;
            for (j = MAX(0, -k); j < tmploop; j++)
                b[i] += a[i][j] * x[j + k];
        }
    }

    public static void cyclic(final double[] a, final double[] b, final double[] c, final double alpha, final double beta, final double[] r,
            final double[] x) throws NRException {
        // Solves for a vector x[0..n-1] the “cyclic” set of linear equations
        // given by equation (2.7.9). a, b, c, and r are input vectors, all
        // dimensioned as [0..n-1], while alpha and beta are the corner entries
        // in the matrix. The input is not modified.
        int i, n = a.length;
        double fact, gamma;
        if (n <= 2)
            throw new NRException("n too small in cyclic");
        final double[] bb = doub_arr(n), u = doub_arr(n), z = doub_arr(n);
        gamma = -b[0]; // Avoid subtraction error in forming bb[0].
        bb[0] = b[0] - gamma; // Set up the diagonal of the modified tridi
        bb[n - 1] = b[n - 1] - alpha * beta / gamma; // agonal system.
        for (i = 1; i < n - 1; i++)
            bb[i] = b[i];
        tridag(a, bb, c, r, x); // Solve A  x D r.
        u[0] = gamma; // Set up the vector u.
        u[n - 1] = alpha;
        for (i = 1; i < n - 1; i++)
            u[i] = 0.0;
        tridag(a, bb, c, u, z); // Solve A  z D u.
        fact = (x[0] + beta * x[n - 1] / gamma) / // Form v  x=.1Cv  z/.
                (1.0 + z[0] + beta * z[n - 1] / gamma);
        for (i = 0; i < n; i++)
            x[i] -= fact * z[i]; // Now get the solution vector x.
    }

    public static void vander(final double[] x, final double[] w, final double[] q) {
        // Solves the Vandermonde linear system PN1
        // iD0 xk i wi D qk .k D 0; : : : ; N  1/. Input consists
        // of the vectors x[0..n-1] and q[0..n-1]; the vector w[0..n-1] is
        // output.
        int i, j, k, n = q.length;
        double b, s, t, xx;
        final double[] c = doub_arr(n);
        if (n == 1)
            w[0] = q[0];
        else {
            for (i = 0; i < n; i++)
                c[i] = 0.0; // Initialize array.
            c[n - 1] = -x[0]; // Coefficients of the master polynomial are found
            for (i = 1; i < n; i++) { // by recursion.
                xx = -x[i];
                for (j = (n - 1 - i); j < (n - 1); j++)
                    c[j] += xx * c[j + 1];
                c[n - 1] += xx;
            }
            for (i = 0; i < n; i++) { // Each subfactor in turn
                xx = x[i];
                t = b = 1.0;
                s = q[n - 1];
                for (k = n - 1; k > 0; k--) { // is synthetically divided,
                    b = c[k] + xx * b;
                    s += q[k - 1] * b; // matrix-multiplied by the right-hand
                                       // side,
                    t = xx * t + b;
                }
                w[i] = s / t; // and supplied with a denominator.
            }
        }
    }

    public static void toeplz(final double[] r, final double[] x, final double[] y) throws NRException {
        // Solves the Toeplitz system
        // PN1
        // jD0 R.N1Cij/xj D yi .i D 0; : : : ; N  1/. The Toeplitz
        // matrix need not be symmetric. y[0..n-1] and r[0..2*n-2] are input
        // arrays; x[0..n-1] is the output array.

        int j, k, m, m1, m2, n1, n = y.length;
        double pp, pt1, pt2, qq, qt1, qt2, sd, sgd, sgn, shn, sxn;
        n1 = n - 1;
        if (r[n1] == 0.0)
            throw new NRException("toeplz-1 singular principal minor");
        x[0] = y[0] / r[n1]; // Initialize for the recursion.
        if (n1 == 0)
            return;
        final double[] g = doub_arr(n1), h = doub_arr(n1);
        g[0] = r[n1 - 1] / r[n1];
        h[0] = r[n1 + 1] / r[n1];
        for (m = 0; m < n; m++) { // Main loop over the recursion.
            m1 = m + 1;
            sxn = -y[m1]; // Compute numerator and denominator for x from eq.
            sd = -r[n1]; // (2.8.19),
            for (j = 0; j < m + 1; j++) {
                sxn += r[n1 + m1 - j] * x[j];
                sd += r[n1 + m1 - j] * g[m - j];
            }
            if (sd == 0.0)
                throw new NRException("toeplz-2 singular principal minor");
            x[m1] = sxn / sd; // whence x.
            for (j = 0; j < m + 1; j++)
                // Eq. (2.8.16).
                x[j] -= x[m1] * g[m - j];
            if (m1 == n1)
                return;
            sgn = -r[n1 - m1 - 1]; // Compute numerator and denominator for G
                                   // and H,
            shn = -r[n1 + m1 + 1]; // eqs. (2.8.24) and (2.8.23),
            sgd = -r[n1];
            for (j = 0; j < m + 1; j++) {
                sgn += r[n1 + j - m1] * g[j];
                shn += r[n1 + m1 - j] * h[j];
                sgd += r[n1 + j - m1] * h[m - j];
            }
            if (sgd == 0.0)
                throw new NRException("toeplz-3 singular principal minor");
            g[m1] = sgn / sgd; // whence G and H.
            h[m1] = shn / sd;
            k = m;
            m2 = (m + 2) >> 1;
            pp = g[m1];
            qq = h[m1];
            for (j = 0; j < m2; j++) {
                pt1 = g[j];
                pt2 = g[k];
                qt1 = h[j];
                qt2 = h[k];
                g[j] = pt1 - pp * qt2;
                g[k] = pt2 - pp * qt1;
                h[j] = qt1 - qq * pt2;
                h[k--] = qt2 - qq * pt1;
            }
        } // Back for another recurrence.
        throw new NRException("toeplz - should not arrive here!");
    }
    
}
