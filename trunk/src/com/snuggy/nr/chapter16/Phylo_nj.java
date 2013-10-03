
package com.snuggy.nr.chapter16;

import static com.snuggy.nr.refs.Refs.*;
import static com.snuggy.nr.util.Static.*;

import com.snuggy.nr.refs.*;

public class Phylo_nj extends Phylagglom {

    // Derived class implementing the neighbor joining (NJ) method.
    private double[] u;

    @Override
    public void premin(final double[][] d, final int[] nextp) {
        // Before finding the minimum we (re-)calculate the u’s.
        int i, j, ncurr = 0;
        double sum;
        for (i = 0; i >= 0; i = nextp[i])
            ncurr++; // Count live entries.
        for (i = 0; i >= 0; i = nextp[i]) { // Compute u[i].
            sum = 0.;
            for (j = 0; j >= 0; j = nextp[j])
                if (i != j)
                    sum += d[i][j];
            u[i] = sum / (ncurr - 2);
        }
    }

    @Override
    public double dminfn(final double[][] d, final int i, final int j) {
        return d[i][j] - u[i] - u[j]; // NJ finds min of this.
    }

    @Override
    public double dbranchfn(final double[][] d, final int i, final int j) {
        return 0.5 * (d[i][j] + u[i] - u[j]); // NJ setting for branch lengths.
    }

    @Override
    public double dnewfn(final double[][] d, final int k, final int i, final int j, final int ni, final int nj) {
        return 0.5 * (d[i][k] + d[j][k] - d[i][j]); // NJ new distances.
    }

    @Override
    public void drootbranchfn(final double[][] d, final int i, final int j, final int ni, final int nj,
            final $double bi, final $double bj) {
        // Since NJ is unrooted, it is a matter of taste how to assign branch
        // lengths to the root.
        // This prescription plots aesthetically.
        $(bi, d[i][j] * (nj - 1 + 1.e-15) / (ni + nj - 2 + 2.e-15));
        $(bj, d[i][j] * (ni - 1 + 1.e-15) / (ni + nj - 2 + 2.e-15));
    }

    Phylo_nj(final double[][] dist) throws InstantiationException, IllegalAccessException {
        this(dist, -1);
    }

    Phylo_nj(final double[][] dist, final int fsr) throws InstantiationException, IllegalAccessException {
        super(dist, fsr);
        u = doub_arr(n);
        makethetree(dist);
    }

}
