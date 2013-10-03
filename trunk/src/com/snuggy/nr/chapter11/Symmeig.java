
package com.snuggy.nr.chapter11;

import static com.snuggy.nr.chapter11.Static.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Symmeig {

    // Computes all eigenvalues and eigenvectors of a real symmetric
    // matrix by reduction to tridiagonal form followed by QL iteration.
    private int n;
    private double[][] z;
    private double[] d, e;
    private boolean yesvecs;

    public double[] d() {
        return d;
    }

    public double[][] z() {
        return z;
    }

    public Symmeig(final double[][] a) throws NRException {
        this(a, true);
    }

    public Symmeig(final double[][] a, final boolean yesvec) throws NRException {
        // Computes all eigenvalues and eigenvectors of a real symmetric matrix
        // a[0..n-1][0..n-1] by reduction to tridiagonal form followed by QL
        // iteration. On output, d[0..n-1] contains the eigenvalues of a sorted
        // into descending order, while z[0..n-1][0..n-1] is a matrix whose
        // columns contain the corresponding normalized eigenvectors. If
        // yesvecs is input as true (the default), then the eigenvectors are
        // computed. If yesvecs is input as false, only the eigenvalues are
        // computed.

        n = (nrows(a));
        z = doub_mat(a);
        d = doub_arr(n);
        e = doub_arr(n);
        yesvecs = (yesvec);

        tred2(); // Reduction to tridiagonal form; see 11.3.
        tqli(); // Eigensystem of tridiagonal matrix; see 11.4.
        sort();
    }

    public Symmeig(final double[] dd, final double[] ee) throws NRException {
        this(dd, ee, true);
    }

    public Symmeig(final double[] dd, final double[] ee, final boolean yesvec) throws NRException {
        // Computes all eigenvalues and (optionally) eigenvectors of a real,
        // symmetric, tridiagonal matrix by QL iteration. On input, dd[0..n-1]
        // contains the diagonal elements of the tridiagonal matrix. The
        // vector ee[0..n-1] inputs the subdiagonal elements of the
        // tridiagonal matrix, with ee[0] arbitrary. Output is the same as
        // the constructor above.

        n = (dd.length);
        d = (dd);
        e = (ee);
        z = doub_mat(n, n);
        yesvecs = (yesvec);

        for (int i = 0; i < n; i++)
            z[i][i] = 1.0;
        tqli();
        sort();
    }

    public void sort() {
        if (yesvecs)
            eigsrt(d, z);
        else
            eigsrt(d);
    }

    public void tred2() {
        // Householder reduction of a real symmetric matrix z[0..n-1][0..n-1].
        // (The input matrix A to Symmeig is stored in z.) On output, z is
        // replaced by the orthogonal matrix Q effecting the transformation.
        // d[0..n-1] contains the diagonal elements of the tridiagonal matrix
        // and e[0..n-1] the off-diagonal elements, with e[0]=0. If yesvecs is
        // false, so that only eigenvalues will subsequently be determined,
        // several statements are omitted, in which case z contains no useful
        // information on output.

        int l, k, j, i;
        double scale, hh, h, g, f;
        for (i = n - 1; i > 0; i--) {
            l = i - 1;
            h = scale = 0.0;
            if (l > 0) {
                for (k = 0; k < i; k++)
                    scale += abs(z[i][k]);
                if (scale == 0.0) // Skip transformation.
                    e[i] = z[i][l];
                else {
                    for (k = 0; k < i; k++) {
                        z[i][k] /= scale; // Use scaled a’s for transformation.
                        h += z[i][k] * z[i][k]; // Form  in h.
                    }
                    f = z[i][l];
                    g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
                    e[i] = scale * g;
                    h -= f * g; // Now h is equation (11.3.4).
                    z[i][l] = f - g; // Store u in row i of z.
                    f = 0.0;
                    for (j = 0; j < i; j++) {
                        if (yesvecs) // Store u=H in column i of z.
                            z[j][i] = z[i][j] / h;
                        g = 0.0; // Form an element of A  u in g.
                        for (k = 0; k < j + 1; k++)
                            g += z[j][k] * z[i][k];
                        for (k = j + 1; k < i; k++)
                            g += z[k][j] * z[i][k];
                        e[j] = g / h; // Form element of p in temporarily unused
                        f += e[j] * z[i][j]; // element of e.
                    }
                    hh = f / (h + h); // Form K, equation (11.3.11).
                    for (j = 0; j < i; j++) { // Form q and store in e
                                              // overwriting p.
                        f = z[i][j];
                        e[j] = g = e[j] - hh * f;
                        for (k = 0; k < j + 1; k++)
                            // Reduce z, equation (11.3.13).
                            z[j][k] -= (f * e[k] + g * z[i][k]);
                    }
                }
            } else
                e[i] = z[i][l];
            d[i] = h;
        }
        if (yesvecs)
            d[0] = 0.0;
        e[0] = 0.0;
        for (i = 0; i < n; i++) { // Begin accumulation of transformation ma
            if (yesvecs) { // trices.
                if (d[i] != 0.0) { // This block skipped when i=0.
                    for (j = 0; j < i; j++) {
                        g = 0.0;
                        for (k = 0; k < i; k++)
                            // Use u and u=H stored in z to form PQ.
                            g += z[i][k] * z[k][j];
                        for (k = 0; k < i; k++)
                            z[k][j] -= g * z[k][i];
                    }
                }
                d[i] = z[i][i];
                z[i][i] = 1.0; // Reset row and column of z to identity
                for (j = 0; j < i; j++)
                    z[j][i] = z[i][j] = 0.0; // matrix for next iteration.
            } else {
                d[i] = z[i][i]; // Only this statement remains.
            }
        }
    }

    public void tqli() throws NRException {
        // QL algorithm with implicit shifts to determine the eigenvalues and
        // (optionally) the eigenvectors of a real, symmetric, tridiagonal
        // matrix, or of a real symmetric matrix previously reduced by
        // tred2 (11.3). On input, d[0..n-1] contains the diagonal elements
        // of the tridiagonal matrix. On output, it returns the eigenvalues.
        // The vector e[0..n-1] inputs the subdiagonal elements of the
        // tridiagonal matrix, with e[0] arbitrary. On output e is destroyed.
        // If the eigenvectors of a tridiagonal matrix are desired, the
        // matrix z[0..n-1][0..n-1] is input as the identity matrix. If the
        // eigenvectors of a matrix that has been reduced by tred2 are required,
        // then z is input as the matrix output by tred2. In either case,
        // column k of z returns the normalized eigenvector corresponding to
        // d[k].

        int m, l, iter, i, k;
        double s, r, p, g, f, dd, c, b;
        final double EPS = EPS(); // numeric_limits<Doub>::epsilon();
        for (i = 1; i < n; i++)
            e[i - 1] = e[i]; // Convenient to renumber the el
        e[n - 1] = 0.0; // ements of e.
        for (l = 0; l < n; l++) {
            iter = 0;
            do {
                for (m = l; m < n - 1; m++) { // Look for a single small
                                              // subdiagonal
                    // element to split the
                    // matrix.
                    dd = abs(d[m]) + abs(d[m + 1]);
                    if (abs(e[m]) <= EPS * dd)
                        break;
                }
                if (m != l) {
                    if (iter++ == 30)
                        throw new NRException("Too many iterations in tqli");
                    g = (d[l + 1] - d[l]) / (2.0 * e[l]); // Form shift.
                    r = pythag(g, 1.0);
                    g = d[m] - d[l] + e[l] / (g + SIGN(r, g)); // This is dm 
                                                               // ks.
                    s = c = 1.0;
                    p = 0.0;
                    for (i = m - 1; i >= l; i--) { // A plane rotation as in the
                                                   // original
                        // QL, followed by Givens
                        // rotations to restore tridiagonal
                        // form.
                        f = s * e[i];
                        b = c * e[i];
                        e[i + 1] = (r = pythag(f, g));
                        if (r == 0.0) { // Recover from underflow.
                            d[i + 1] -= p;
                            e[m] = 0.0;
                            break;
                        }
                        s = f / r;
                        c = g / r;
                        g = d[i + 1] - p;
                        r = (d[i] - g) * s + 2.0 * c * b;
                        d[i + 1] = g + (p = s * r);
                        g = c * r - b;
                        if (yesvecs) {
                            for (k = 0; k < n; k++) { // Form eigenvectors.
                                f = z[k][i + 1];
                                z[k][i + 1] = s * z[k][i] + c * f;
                                z[k][i] = c * z[k][i] - s * f;
                            }
                        }
                    }
                    if (r == 0.0 && i >= l)
                        continue;
                    d[l] -= p;
                    e[l] = g;
                    e[m] = 0.0;
                }
            } while (m != l);
        }
    }

    public double pythag(final double a, final double b) {
        // Computes .a2 Cb2/1=2 without destructive underflow or overflow.
        double absa = abs(a), absb = abs(b);
        return (absa > absb ? absa * sqrt(1.0 + SQR(absb / absa)) : (absb == 0.0 ? 0.0 : absb
                * sqrt(1.0 + SQR(absa / absb))));
    }

    // void tred2();
    // void tqli();
    // Doub pythag(final double a, final double b);

}
