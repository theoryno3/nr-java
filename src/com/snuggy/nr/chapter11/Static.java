package com.snuggy.nr.chapter11;

public class Static {

    public static void eigsrt(final double[] d) {
        eigsrt(d, null);
    }

    public static void eigsrt(final double[] d, final double[][] v) {
        // Given the eigenvalues d[0..n-1] and (optionally) the eigenvectors
        // v[0..n-1][0..n-1] as determined by Jacobi (11.1) or tqli (11.4),
        // this routine sorts the eigenvalues into descending order and
        // rearranges the columns of v correspondingly. The method is
        // straight insertion.

        int k;
        int n = d.length;
        for (int i = 0; i < n - 1; i++) {
            double p = d[k = i];
            for (int j = i; j < n; j++)
                if (d[j] >= p)
                    p = d[k = j];
            if (k != i) {
                d[k] = d[i];
                d[i] = p;
                if (v != null)
                    for (int j = 0; j < n; j++) {
                        p = v[j][i];
                        v[j][i] = v[j][k];
                        v[j][k] = p;
                    }
            }
        }
    }
    
}
