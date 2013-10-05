package com.snuggy.nr.chapter02;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class SVD {

    // Object for singular value decomposition of a matrix A, and related
    // functions.
    private int m, n;
    private final double[][] u, v; // The matrices U and V.
    private final double[] w; // The diagonal matrix W.
    private double eps, tsh;

    public SVD(final double[][] a) throws NRException {
        // Constructor. The single argument is A. The SVD computation is done
        // by decompose, and the results are sorted by reorder.
        m = (nrows(a));
        n = (ncols(a));
        u = doub_mat(a);
        v = doub_mat(n, n);
        w = doub_vec(n);
        eps = EPS();
        decompose();
        reorder();
        tsh = 0.5 * sqrt(m + n + 1.) * w[0] * eps; // Default threshold for
                                                        // nonzero singular
    } // values.
    
    public double tsh() {
        return tsh;
    }
    
    public final double[][] v() {
        return v;
    }
    
    public final double[] w() {
        return w;
    }

    // void solve(VecDoub_I &b, VecDoub_O &x, Doub thresh);
    // void solve(MatDoub_I &b, MatDoub_O &x, Doub thresh);

    // Solve with (apply the pseudoinverse to) one or more right-hand sides.
    // Int rank(Doub thresh); Quantities associated with the range and
    // Int nullity(Doub thresh); nullspace of A.
    // MatDoub range(Doub thresh);
    // MatDoub nullspace(Doub thresh);

    public double inv_condition() { // Return reciprocal of the condition num
        return (w[0] <= 0. || w[n - 1] <= 0.) ? 0. : w[n - 1] / w[0]; // ber of
                                                                      // A.
    }

    public void decompose() throws NRException {
        boolean flag;
        int i = 0;
        int its = 0;
        int j = 0;
        int jj = 0;
        int k = 0;
        int l = 0;
        int nm = 0;
        double anorm = 0.0;
        double c = 0.0;
        double f = 0.0;
        double g = 0.0;
        double h = 0.0;
        double s = 0.0;
        double scale = 0.0;
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        final double[] rv1 = doub_vec(n);
        g = scale = anorm = 0.0;
        for (i = 0; i < n; i++) {
            l = i + 2;
            rv1[i] = scale * g;
            g = s = scale = 0.0;
            if (i < m) {
                for (k = i; k < m; k++)
                    scale += abs(u[k][i]);
                if (scale != 0.0) {
                    for (k = i; k < m; k++) {
                        u[k][i] /= scale;
                        s += u[k][i] * u[k][i];
                    }
                    f = u[i][i];
                    g = -SIGN(sqrt(s), f);
                    h = f * g - s;
                    u[i][i] = f - g;
                    for (j = l - 1; j < n; j++) {
                        for (s = 0.0, k = i; k < m; k++)
                            s += u[k][i] * u[k][j];
                        f = s / h;
                        for (k = i; k < m; k++)
                            u[k][j] += f * u[k][i];
                    }
                    for (k = i; k < m; k++)
                        u[k][i] *= scale;
                }
            }
            w[i] = scale * g;
            g = s = scale = 0.0;
            if (i + 1 <= m && i + 1 != n) {
                for (k = l - 1; k < n; k++)
                    scale += abs(u[i][k]);
                if (scale != 0.0) {
                    for (k = l - 1; k < n; k++) {
                        u[i][k] /= scale;
                        s += u[i][k] * u[i][k];
                    }
                    f = u[i][l - 1];
                    g = -SIGN(sqrt(s), f);
                    h = f * g - s;
                    u[i][l - 1] = f - g;
                    for (k = l - 1; k < n; k++)
                        rv1[k] = u[i][k] / h;
                    for (j = l - 1; j < m; j++) {
                        for (s = 0.0, k = l - 1; k < n; k++)
                            s += u[j][k] * u[i][k];
                        for (k = l - 1; k < n; k++)
                            u[j][k] += s * rv1[k];
                    }
                    for (k = l - 1; k < n; k++)
                        u[i][k] *= scale;
                }
            }
            anorm = MAX(anorm, (abs(w[i]) + abs(rv1[i])));
        }
        for (i = n - 1; i >= 0; i--) {
            if (i < n - 1) {
                if (g != 0.0) {
                    for (j = l; j < n; j++)
                        v[j][i] = (u[i][j] / u[i][l]) / g;
                    for (j = l; j < n; j++) {
                        for (s = 0.0, k = l; k < n; k++)
                            s += u[i][k] * v[k][j];
                        for (k = l; k < n; k++)
                            v[k][j] += s * v[k][i];
                    }
                }
                for (j = l; j < n; j++)
                    v[i][j] = v[j][i] = 0.0;
            }
            v[i][i] = 1.0;
            g = rv1[i];
            l = i;
        }
        for (i = MIN(m, n) - 1; i >= 0; i--) {
            l = i + 1;
            g = w[i];
            for (j = l; j < n; j++)
                u[i][j] = 0.0;
            if (g != 0.0) {
                g = 1.0 / g;
                for (j = l; j < n; j++) {
                    for (s = 0.0, k = l; k < m; k++)
                        s += u[k][i] * u[k][j];
                    f = (s / u[i][i]) * g;
                    for (k = i; k < m; k++)
                        u[k][j] += f * u[k][i];
                }
                for (j = i; j < m; j++)
                    u[j][i] *= g;
            } else
                for (j = i; j < m; j++)
                    u[j][i] = 0.0;
            ++u[i][i];
        }
        for (k = n - 1; k >= 0; k--) {
            for (its = 0; its < 30; its++) {
                flag = true;
                for (l = k; l >= 0; l--) {
                    nm = l - 1;
                    if (l == 0 || abs(rv1[l]) <= eps * anorm) {
                        flag = false;
                        break;
                    }
                    if (abs(w[nm]) <= eps * anorm)
                        break;
                }
                if (flag) {
                    c = 0.0;
                    s = 1.0;
                    for (i = l; i < k + 1; i++) {
                        f = s * rv1[i];
                        rv1[i] = c * rv1[i];
                        if (abs(f) <= eps * anorm)
                            break;
                        g = w[i];
                        h = pythag(f, g);
                        w[i] = h;
                        h = 1.0 / h;
                        c = g * h;
                        s = -f * h;
                        for (j = 0; j < m; j++) {
                            y = u[j][nm];
                            z = u[j][i];
                            u[j][nm] = y * c + z * s;
                            u[j][i] = z * c - y * s;
                        }
                    }
                }
                z = w[k];
                if (l == k) {
                    if (z < 0.0) {
                        w[k] = -z;
                        for (j = 0; j < n; j++)
                            v[j][k] = -v[j][k];
                    }
                    break;
                }
                if (its == 29)
                    throw new NRException("no convergence in 30 svdcmp iterations");
                x = w[l];
                nm = k - 1;
                y = w[nm];
                g = rv1[nm];
                h = rv1[k];
                f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
                g = pythag(f, 1.0);
                f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
                c = s = 1.0;
                for (j = l; j <= nm; j++) {
                    i = j + 1;
                    g = rv1[i];
                    y = w[i];
                    h = s * g;
                    g = c * g;
                    z = pythag(f, h);
                    rv1[j] = z;
                    c = f / z;
                    s = h / z;
                    f = x * c + g * s;
                    g = g * c - x * s;
                    h = y * s;
                    y *= c;
                    for (jj = 0; jj < n; jj++) {
                        x = v[jj][j];
                        z = v[jj][i];
                        v[jj][j] = x * c + z * s;
                        v[jj][i] = z * c - x * s;
                    }
                    z = pythag(f, h);
                    w[j] = z;
                    if (z != 0.0) {
                        z = 1.0 / z;
                        c = f * z;
                        s = h * z;
                    }
                    f = c * g + s * y;
                    x = c * y - s * g;
                    for (jj = 0; jj < m; jj++) {
                        y = u[jj][j];
                        z = u[jj][i];
                        u[jj][j] = y * c + z * s;
                        u[jj][i] = z * c - y * s;
                    }
                }
                rv1[l] = 0.0;
                rv1[k] = f;
                w[k] = x;
            }
        }
    }

    public void reorder() {
        int i = 0;
        int j = 0;
        int k = 0;
        int s = 0;
        int inc = 1;
        double sw = 0.0;
        final double[] su = doub_vec(m);
        final double[] sv = doub_vec(n);
        do {
            inc *= 3;
            inc++;
        } while (inc <= n);
        do {
            inc /= 3;
            for (i = inc; i < n; i++) {
                sw = w[i];
                for (k = 0; k < m; k++)
                    su[k] = u[k][i];
                for (k = 0; k < n; k++)
                    sv[k] = v[k][i];
                j = i;
                while (w[j - inc] < sw) {
                    w[j] = w[j - inc];
                    for (k = 0; k < m; k++)
                        u[k][j] = u[k][j - inc];
                    for (k = 0; k < n; k++)
                        v[k][j] = v[k][j - inc];
                    j -= inc;
                    if (j < inc)
                        break;
                }
                w[j] = sw;
                for (k = 0; k < m; k++)
                    u[k][j] = su[k];
                for (k = 0; k < n; k++)
                    v[k][j] = sv[k];

            }
        } while (inc > 1);
        for (k = 0; k < n; k++) {
            s = 0;
            for (i = 0; i < m; i++)
                if (u[i][k] < 0.)
                    s++;
            for (j = 0; j < n; j++)
                if (v[j][k] < 0.)
                    s++;
            if (s > (m + n) / 2) {
                for (i = 0; i < m; i++)
                    u[i][k] = -u[i][k];
                for (j = 0; j < n; j++)
                    v[j][k] = -v[j][k];
            }
        }
    }

    public double pythag(final double a, final double b) {
        double absa = abs(a);
        double absb = abs(b);
        return (absa > absb ? 
                absa * sqrt(1.0 + SQR(absb / absa)) : 
                (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + SQR(absa / absb))));
    }

    // void decompose(); Functions used by the constructor.
    // void reorder();
    // Doub pythag(final double a, final double b);

    public int rank() {
        return rank(-1.);
    }

    public int rank(final double thresh) {
        // Return the rank of A, after zeroing any singular values smaller than
        // thresh. If thresh is negative, a default value based on estimated
        // roundoff is used.
        int j, nr = 0;
        tsh = (thresh >= 0. ? thresh : 0.5 * sqrt(m + n + 1.) * w[0] * eps);
        for (j = 0; j < n; j++)
            if (w[j] > tsh)
                nr++;
        return nr;
    }

    public int nullity() {
        return nullity(-1.);
    }

    public int nullity(final double thresh) {
        // Return the nullity of A, after zeroing any singular values smaller
        // than thresh. Default value as above.
        int j, nn = 0;
        tsh = (thresh >= 0. ? thresh : 0.5 * sqrt(m + n + 1.) * w[0] * eps);
        for (j = 0; j < n; j++)
            if (w[j] <= tsh)
                nn++;
        return nn;
    }

    public final double[][] range() {
        return range(-1.);
    }

    public final double[][] range(final double thresh) {
        // Give an orthonormal basis for the range of A as the columns of a
        // returned matrix. thresh as above.
        int i, j, nr = 0;
        final double[][] rnge = doub_mat(m, rank(thresh));
        for (j = 0; j < n; j++) {
            if (w[j] > tsh) {
                for (i = 0; i < m; i++)
                    rnge[i][nr] = u[i][j];
                nr++;
            }
        }
        return rnge;
    }

    public final double[][] nullspace() {
        return nullspace(-1.);
    }

    public final double[][] nullspace(final double thresh) {
        // Give an orthonormal basis for the nullspace of A as the columns of
        // a returned matrix. thresh as above.
        int j, jj, nn = 0;
        final double[][] nullsp = doub_mat(n, nullity(thresh));
        for (j = 0; j < n; j++) {
            if (w[j] <= tsh) {
                for (jj = 0; jj < n; jj++)
                    nullsp[jj][nn] = v[jj][j];
                nn++;
            }
        }
        return nullsp;
    }

    public void solve(final double[] b, final double[] x) throws NRException {
        solve(b, x, -1.);
    }

    public void solve(final double[] b, final double[] x, final double thresh) throws NRException {
        // Solve A  x D b for a vector x using the pseudoinverse of A as
        // obtained by SVD. If positive, thresh is the threshold value below
        // which singular values are considered as zero. If thresh is
        // negative, a default based on expected roundoff error is used.
        int i, j, jj;
        double s;
        if (b.length != m || x.length != n)
            throw new NRException("SVD::solve bad sizes");
        final double[] tmp = doub_vec(n);
        tsh = (thresh >= 0. ? thresh : 0.5 * sqrt(m + n + 1.) * w[0] * eps);
        for (j = 0; j < n; j++) { // Calculate UT B.
            s = 0.0;
            if (w[j] > tsh) { // Nonzero result only if wj is nonzero.
                for (i = 0; i < m; i++)
                    s += u[i][j] * b[i];
                s /= w[j]; // This is the divide by wj .
            }
            tmp[j] = s;
        }
        for (j = 0; j < n; j++) { // Matrix multiply by V to get answer.
            s = 0.0;
            for (jj = 0; jj < n; jj++)
                s += v[j][jj] * tmp[jj];
            x[j] = s;
        }
    }

    public void solve(final double[][] b, final double[][] x) throws NRException {
        solve(b, x, -1.);
    }

    public void solve(final double[][] b, final double[][] x, final double thresh) throws NRException {
        // Solves m sets of n equations AX D B using the pseudoinverse of A.
        // The right-hand sides are input as b[0..n-1][0..m-1], while
        // x[0..n-1][0..m-1] returns the solutions. thresh as above.
        int i, j, m = ncols(b);
        if (nrows(b) != n || nrows(x) != n || ncols(b) != ncols(x))
            throw new NRException("SVD::solve bad sizes");
        final double[] xx = doub_vec(n);
        for (j = 0; j < m; j++) { // Copy and solve each column in turn.
            for (i = 0; i < n; i++)
                xx[i] = b[i][j];
            solve(xx, xx, thresh);
            for (i = 0; i < n; i++)
                x[i][j] = xx[i];
        }
    }

}
