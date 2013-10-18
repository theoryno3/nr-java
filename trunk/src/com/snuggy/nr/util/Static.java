
package com.snuggy.nr.util;

import static java.lang.Math.*;

import java.lang.reflect.*;

import com.snuggy.nr.refs.*;

public class Static {

    // VecDoub v --> final double[] v (v[i])
    // MatDoub m --> final double[][] m (m[i][j])
    // VecInt v --> final int[] v (v[i])
    // int& i --> final int[] i (i[0])
    // double& d --> final double[] d (d[0])
    // double* p --> final double[] p_arr, int p_off (p_arr[p_off])

    public static final double[][] mat_ident(int n) {
        final double[][] r;
        r = new double[n][n];
        for (int i = 0; i < n; i++)
            r[i][i] = 1.0;
        return r;
    }

    public static final int[] times(final int[][] a, final int[] b) {
        int m = a.length;
        final int[] r = new int[m];
        for (int i = 0; i < m; i++) {
            int sum = 0;
            for (int k = 0; k < m; k++)
                sum += a[i][k] * b[k];
            r[i] = sum;
        }
        return r;
    }

    public static final double[] times(final double[][] a, final double[] b) {
        int m = a.length;
        final double[] r = new double[m];
        for (int i = 0; i < m; i++) {
            double sum = 0;
            for (int k = 0; k < m; k++)
                sum += a[i][k] * b[k];
            r[i] = sum;
        }
        return r;
    }

    public static final double[][] mat_times(final double[][] a, final double[][] b) throws NRException {
        int m = a.length;
        int n = b[0].length;
        int l = a[0].length;
        if (l != b.length)
            throw new NRException("l != b.length");
        final double[][] r = new double[m][n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++) {
                double sum = 0.0;
                for (int k = 0; k < l; k++)
                    sum += a[i][k] * b[k][j];
                r[i][j] = sum;
            }
        return r;
    }

    public static final double[][] transpose(final double[][] a) {
        int m = a.length;
        int n = a[0].length;
        final double[][] r = new double[n][m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                r[j][i] = a[i][j];
        return r;
    }
    
    @SuppressWarnings("unchecked")
    public static <T> T[] obj_vec(Class<T> type, int n) throws InstantiationException, IllegalAccessException {
        T[] arr = (T[]) Array.newInstance(type, n);
        for (int i = 0; i < n; i++)
            arr[i] = type.newInstance();
        return arr;
    }

    @SuppressWarnings("unchecked")
    public static <T> T[] obj_vec(Class<T> type, final T[] a) throws InstantiationException, IllegalAccessException {
        T[] arr = (T[]) Array.newInstance(type, a.length);
        for (int i = 0; i < a.length; i++)
            arr[i] = a[i];
        return arr;
    }

    public static <T> T[] obj_vec_nulls(Class<T> type, int n) throws InstantiationException, IllegalAccessException {
        @SuppressWarnings("unchecked")
        T[] arr = (T[]) Array.newInstance(type, n);
        for (int i = 0; i < n; i++)
            arr[i] = null;
        return arr;
    }

    public static boolean[] bool_vec(int n) {
        return new boolean[n];
    }

    public static final int[] int_vec(int n) {
        return new int[n];
    }

    public static final int[][] int_mat(int m, int n) {
        return new int[m][n];
    }

    public static final double[] doub_vec(int n) {
        return new double[n];
    }

    public static final double[][] doub_mat(int n, int m) {
        return new double[n][m];
    }

    public static final double[][] doub_mat(int n, int m, double x) {
        final double[][] r = new double[n][m];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                r[i][j] = x;
        return r;
    }

    public static final double[][] doub_mat(int n, int m, final double[] x) throws NRException {
        if (x.length != n * m)
            throw new NRException("x.length != n * m");
        final double[][] r = new double[n][m];
        int p = 0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                r[i][j] = x[p++];
        return r;
    }

    public static final double[] doub_vec(int n, double x) {
        final double[] r = new double[n];
        for (int i = 0; i < n; i++)
            r[i] = x;
        return r;
    }

    public static long[] long_vec(int n) {
        return new long[n];
    }

    public static double frexp(double x, $int exp_out) {
        double y = abs(x);
        int exp = 0;
        while (y < 0.5) {
            y *= 2.0;
            exp--;
        }
        while (y >= 1.0) {
            y /= 2.0;
            exp++;
        }
        exp_out.$(exp);
        return signum(x) * y;
    }

    public static double ldexp(double x, int exp) {
        return x * pow(2.0, exp);
    }

    public static double Doub(int x) {
        return (double) x;
    }

    public static int Int(double x) {
        return (int) x;
    }

    public static int nrows(final double[][] m) {
        return m.length;
    }

    public static int ncols(final double[][] m) {
        return m[0].length;
    }

    public static void SWAP(final double[] x, final int i, final int j) {
        double t = x[i];
        x[i] = x[j];
        x[j] = t;
    }

    public static void SWAP(final int[] x, final int i, final int j) {
        int t = x[i];
        x[i] = x[j];
        x[j] = t;
    }

    public static void SWAP(final $$int1d x, final int i, final int j) {
        int t = x.$()[i];
        x.$()[i] = x.$()[j];
        x.$()[j] = t;
    }

    public static <T> void SWAP(final T[] x, final int i, final int j) {
        T t = x[i];
        x[i] = x[j];
        x[j] = t;
    }

    public static void SWAP(final double[][] x, final int i, final int j, final int k, final int l) {
        double t = x[i][j];
        x[i][j] = x[k][l];
        x[k][l] = t;
    }

    public static int MIN(final int a, final int b) {
        return (a < b ? a : b);
    }

    public static int MAX(final int a, final int b) {
        return (a > b ? a : b);
    }

    public static double MIN(final double a, final double b) {
        return (a < b ? a : b);
    }

    public static double MAX(final double a, final double b) {
        return (a > b ? a : b);
    }

    public static double SQR(final double a) {
        return a * a;
    }

    public static int SQR(final int a) {
        return a * a;
    }

    public static final double[][] doub_mat(final double[][] m) {
        final double[][] r = new double[m.length][];
        for (int i = 0; i < m.length; i++) {
            r[i] = new double[m[i].length];
            System.arraycopy(m[i], 0, r[i], 0, m[i].length);
        }
        return r;
    }

    public static final int[][] int_mat(final int[][] m) {
        final int[][] r = new int[m.length][];
        for (int i = 0; i < m.length; i++) {
            r[i] = new int[m[i].length];
            System.arraycopy(m[i], 0, r[i], 0, m[i].length);
        }
        return r;
    }

    public static final int[] int_vec(final int n, final int x) {
        final int[] r = new int[n];
        for (int i = 0; i < n; i++)
            r[i] = x;
        return r;
    }

    public static final int[] int_vec(final int n, final int[] from, final int offset) {
        final int[] r = new int[n];
        for (int i = 0; i < n; i++)
            r[i] = from[offset + i];
        return r;
    }

    public static final int[][] int_mat(final int m, final int n, final int[] a) {
        final int[][] r = new int[m][n];
        int p = 0;
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                r[i][j] = a[p++];
        return r;
    }

    public static final double[] doub_vec(final double[] a) {
        final double[] r = new double[a.length];
        System.arraycopy(a, 0, r, 0, a.length);
        return r;
    }

    public static final int[] int_vec(final int[] a) {
        final int[] r = new int[a.length];
        System.arraycopy(a, 0, r, 0, a.length);
        return r;
    }

    public static void main(String[] argv) {
        double x = 1.1;
        while (x > 1.0) {
            System.out.println(x - 1.0);
            x = (x - 1.0) / 3.0 + 1.0;
        }
    }

    public static double EPS() {
        return 2.220446049250313E-16;
    }

    public static double SIGN(final double a, final double b) {
        return (b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));
    }

    public static double max_norm(final double[] arr) {
        double r = 0.0;
        for (double x : arr)
            if (abs(x) > r)
                r = abs(x);
        return r;
    }

    public static int max_norm(final int[] arr) {
        int r = 0;
        for (int x : arr)
            if (abs(x) > r)
                r = abs(x);
        return r;
    }

    public static double max_norm(final double[][] arr) {
        double r = 0.0;
        for (final double[] row : arr)
            for (double x : row)
                if (abs(x) > r)
                    r = abs(x);
        return r;
    }

    public static final double[] vec_minus(final double[] a, final double[] b) throws NRException {
        if (a == null || b == null || a.length != b.length)
            throw new NRException("a == null || b == null || a.length != b.length");
        final double[] r = new double[a.length];
        for (int i = 0; i < a.length; i++)
            r[i] = a[i] - b[i];
        return r;
    }

    public static final int[] vec_minus(final int[] a, final int[] b) throws NRException {
        if (a == null || b == null || a.length != b.length)
            throw new NRException("a == null || b == null || a.length != b.length");
        final int[] r = new int[a.length];
        for (int i = 0; i < a.length; i++)
            r[i] = a[i] - b[i];
        return r;
    }

    public static final double[] vec_plus(final double[] a, final double[] b) throws NRException {
        if (a == null || b == null || a.length != b.length)
            throw new NRException("a == null || b == null || a.length != b.length");
        final double[] r = new double[a.length];
        for (int i = 0; i < a.length; i++)
            r[i] = a[i] + b[i];
        return r;
    }

    public static final double[][] mat_minus(final double[][] a, final double[][] b) throws NRException {
        int i, j, m = a.length, n = a[0].length;
        if (a.length != b.length || a[0].length != b[0].length)
            throw new NRException("a.length != b.length || a[0].length != b[0].length");
        final double[][] c = new double[m][n];
        for (i = 0; i < m; i++)
            for (j = 0; j < n; j++) {
                c[i][j] = a[i][j] - b[i][j];
            }
        return c;
    }

    public static final double[][] mat_plus(final double[][] a, final double[][] b) throws NRException {
        int i, j, m = a.length, n = a[0].length;
        if (a.length != b.length || a[0].length != b[0].length)
            throw new NRException("a.length != b.length || a[0].length != b[0].length");
        final double[][] c = new double[m][n];
        for (i = 0; i < m; i++)
            for (j = 0; j < n; j++) {
                c[i][j] = a[i][j] + b[i][j];
            }
        return c;
    }
    
    public static long umod(long x, int m) throws NRException {
        long r = ((x % m) + m) % m;
        return r;
    }
    
    public static int umod(int x, int m) throws NRException {
        int r = ((x % m) + m) % m;
        return r;
    }
    
    /*
     * private static double cof[] = { 76.18009173, -86.50532033, 24.01409822,
     * -1.231739516, 0.120858003e-2, -0.536382e-5 };
     * 
     * public static double gammln(final double xx) { double x, tmp, ser; int j;
     * 
     * x = xx - 1.0; tmp = x + 5.5; tmp -= (x + 0.5) * log(tmp); ser = 1.0; for
     * (j = 0; j <= 5; j++) { x += 1.0; ser += cof[j] / x; } return -tmp +
     * log(2.50662827465 * ser); }
     */
}
