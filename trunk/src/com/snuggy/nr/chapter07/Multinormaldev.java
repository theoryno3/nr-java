
package com.snuggy.nr.chapter07;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.chapter02.*;
import com.snuggy.nr.util.*;

public class Multinormaldev extends Ran {

    // Structure for multivariate normal deviates.
    private int mm;
    private double[] mean;
    private double[][] var;
    private Cholesky chol;
    private double[] spt, pt;

    public Multinormaldev(final long j, final double[] mmean, final double[][] vvar) throws NRException {
        // Constructor. Arguments are the random generator seed, the (vector)
        // mean, and the (matrix) covariance. Cholesky decomposition of the
        // covariance is done here.
        super(j);
        mm = (mmean.length);
        mean = (mmean);
        var = (vvar);
        chol = new Cholesky(var);
        spt = doub_arr(mm);
        pt = doub_arr(mm);
        if (ncols(var) != mm || nrows(var) != mm)
            throw new NRException("bad sizes");
    }

    public double[] dev() throws NRException {
        // Return a multivariate normal deviate.
        int i;
        double u, v, x, y, q;
        for (i = 0; i < mm; i++) { // Fill a vector of independent normal
                                   // deviates.
            do {
                u = doub();
                v = 1.7156 * (doub() - 0.5);
                x = u - 0.449871;
                y = abs(v) + 0.386595;
                q = SQR(x) + y * (0.19600 * y - 0.25472 * x);
            } while (q > 0.27597 && (q > 0.27846 || SQR(v) > -4. * log(u) * SQR(u)));
            spt[i] = v / u;
        }
        chol.elmult(spt, pt); // Apply equation (7.4.3).
        for (i = 0; i < mm; i++) {
            pt[i] += mean[i];
        }
        return pt;
    }
}
