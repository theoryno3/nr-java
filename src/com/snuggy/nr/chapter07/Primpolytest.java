
package com.snuggy.nr.chapter07;

import static com.snuggy.nr.util.Static.*;

import com.snuggy.nr.util.*;

public class Primpolytest {

    // Test polynomials over the integers mod 2 for primitiveness.
    private int N, nfactors;
    private long[] factors;
    private final int[] t, a, p;

    public Primpolytest() {
        // Constructor. The constants are specific to 32-bit LFSRs.
        N = (32);
        nfactors = (5);
        factors = long_vec(nfactors);
        t = int_vec(N * N);
        a = int_vec(N * N);
        p = int_vec(N * N);
        long factordata[] = { 3, 5, 17, 257, 65537 };
        for (int i = 0; i < nfactors; i++)
            factors[i] = factordata[i];
    }

    public int ispident() {
        // Utility to test if p is the identity matrix.
        int i, j;
        for (i = 0; i < N; i++)
            for (j = 0; j < N; j++) {
                if (i == j) {
                    if (p[i * N + j] != 1)
                        return 0;
                } else {
                    if (p[i * N + j] != 0)
                        return 0;
                }
            }
        return 1;
    }

    public void mattimeseq(final int[] a, final int[] b) {
        // Utility for a *= b on matrices a and b.
        int i, j, k, sum;
        final int[] tmp = int_vec(N * N);
        for (i = 0; i < N; i++)
            for (j = 0; j < N; j++) {
                sum = 0;
                for (k = 0; k < N; k++)
                    sum += a[i * N + k] * b[k * N + j];
                tmp[i * N + j] = sum & 1;
            }
        for (k = 0; k < N * N; k++)
            a[k] = tmp[k];
    }

    public void matpow(long n) {
        // Utility for matrix p = a^n by successive
        int k; // squares.
        for (k = 0; k < N * N; k++)
            p[k] = 0;
        for (k = 0; k < N; k++)
            p[k * N + k] = 1;
        while (true) {
            if ((n & 1) != 0)
                mattimeseq(p, a);
            n >>= 1;
            if (n == 0)
                break;
            mattimeseq(a, a);
        }
    }

    public int test(long n) throws NRException {
        // Main test routine. Returns 1 if the polynomial with serial number n
        // (see text) is primitive, 0 otherwise.
        int i, k, j;
        long pow, tnm1, nn = n;
        tnm1 = ((long) 1 << N) - 1;
        if (n > (tnm1 >>> 1))
            throw new NRException("not a polynomial of degree N");
        for (k = 0; k < N * N; k++)
            t[k] = 0; // Construct the update matrix in t.
        for (i = 1; i < N; i++)
            t[i * N + (i - 1)] = 1;
        j = 0;
        while (nn != 0) {
            if ((nn & 1) != 0)
                t[j] = 1;
            nn >>= 1;
            j++;
        }
        t[N - 1] = 1;
        for (k = 0; k < N * N; k++)
            a[k] = t[k]; // Test that t^tnm1 is the identity matrix.
        matpow(tnm1);
        if (ispident() != 1)
            return 0;
        for (i = 0; i < nfactors; i++) { // Test that the t to the required
                                         // submulti
            pow = tnm1 / factors[i]; // ple powers is not the identity matrix.
            for (k = 0; k < N * N; k++)
                a[k] = t[k];
            matpow(pow);
            if (ispident() == 1)
                return 0;
        }
        return 1;
    }
}
