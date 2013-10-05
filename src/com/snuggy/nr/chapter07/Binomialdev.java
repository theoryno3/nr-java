package com.snuggy.nr.chapter07;

import static com.snuggy.nr.chapter06.Static.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Binomialdev extends Ran {

    // Structure for binomial deviates.
    private double pp, p, pb, np, glnp, plog, pclog, sq;
    @SuppressWarnings("unused")
    private double expnp;
    private int n, swch;
    private long uz, uo, unfin, diff, rltp;
    private final int[] pbits = int_arr(5);
    private final double[] cdf = doub_arr(64);
    private final double[] logfact = doub_arr(1024);

    public Binomialdev(final int nn, final double ppp, final long i) throws NRException {
        // Constructor arguments are n, p, and a random sequence seed.
        super(i);
        pp = (ppp);
        n = (nn);
        int j;
        pb = p = (pp <= 0.5 ? pp : 1.0 - pp);
        if (n <= 64) { // Will use bit-parallel direct method.
            uz = 0;
            uo = 0xffffffffffffffffL;
            rltp = 0;
            for (j = 0; j < 5; j++)
                pbits[j] = 1 & ((int) (pb *= 2.));
            pb -= floor(pb); // Leading bits of p (above) and remaining
            swch = 0; // fraction.
        } else if (n * p < 30.) { // Will use precomputed cdf table.
            cdf[0] = exp(n * log(1 - p));
            for (j = 1; j < 64; j++)
                cdf[j] = cdf[j - 1]
                        + exp(gammln(n + 1.) - gammln(j + 1.) - gammln(n - j + 1.) + j * log(p) + (n - j) * log(1. - p));
            swch = 1;
        } else { // Will use ratio-of-uniforms method.
            np = n * p;
            glnp = gammln(n + 1.);
            plog = log(p);
            pclog = log(1. - p);
            sq = sqrt(np * (1. - p));
            if (n < 1024)
                for (j = 0; j <= n; j++)
                    logfact[j] = gammln(j + 1.);
            swch = 2;
        }
    }

    public int dev() throws NRException {
        // Return a binomial deviate.
        int j, k, kl, km;
        double y, u, v, u2, v2, b;
        if (swch == 0) {
            unfin = uo; // Mark all bits as ”unfinished.”
            for (j = 0; j < 5; j++) { // Compare with first five bits of p.
                diff = unfin & (int64() ^ ((pbits[j] != 0) ? uo : uz)); // Mask
                                                                        // of
                // diff.
                if (pbits[j] != 0)
                    rltp |= diff; // Set bits to 1, meaning ran < p.
                else
                    rltp = rltp & ~diff; // Set bits to 0, meaning ran > p.
                unfin = unfin & ~diff; // Update unfinished status.
            }
            k = 0; // Now we just count the events.
            for (j = 0; j < n; j++) {
                if ((unfin & 1) != 0) {
                    if (doub() < pb)
                        ++k;
                } // Clean up unresolved cases,
                else {
                    if ((rltp & 1) != 0)
                        ++k;
                } // or use bit answer.
                unfin >>= 1;
                rltp >>= 1;
            }
        } else if (swch == 1) { // Use stored cdf.
            y = doub();
            kl = -1;
            k = 64;
            while (k - kl > 1) {
                km = (kl + k) / 2;
                if (y < cdf[km])
                    k = km;
                else
                    kl = km;
            }
        } else { // Use ratio-of-uniforms method.
            for (;;) {
                u = 0.645 * doub();
                v = -0.63 + 1.25 * doub();
                v2 = SQR(v);
                // Try squeeze for fast rejection:
                if (v >= 0.) {
                    if (v2 > 6.5 * u * (0.645 - u) * (u + 0.2))
                        continue;
                } else {
                    if (v2 > 8.4 * u * (0.645 - u) * (u + 0.1))
                        continue;
                }
                k = Int(floor(sq * (v / u) + np + 0.5));
                if (k < 0)
                    continue;
                u2 = SQR(u);
                // Try squeeze for fast acceptance:
                if (v >= 0.) {
                    if (v2 < 12.25 * u2 * (0.615 - u) * (0.92 - u))
                        break;
                } else {
                    if (v2 < 7.84 * u2 * (0.615 - u) * (1.2 - u))
                        break;
                }
                b = sq * exp(glnp + k * plog + (n - k) * pclog // Only when we
                                                               // must.
                        - (n < 1024 ? logfact[k] + logfact[n - k] : gammln(k + 1.) + gammln(n - k + 1.)));
                if (u2 < b)
                    break;
            }
        }
        if (p != pp)
            k = n - k;
        return k;
    }

}
