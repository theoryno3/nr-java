
package com.snuggy.nr.chapter16;
import static com.snuggy.nr.refs.Refs.*;

import com.snuggy.nr.refs.*;

public class Phylo_upgma extends Phylagglom {

    // Derived class implementing the UPGMA method. Only need to define
    // functions that are virtual in Phylagglom.

    @Override
    public void premin(final double[][] d, final int[] nextp) {
        // No pre-min calculations.
    }

    @Override
    public double dminfn(final double[][] d, final int i, final int j) {
        return d[i][j];
    }

    @Override
    public double dbranchfn(final double[][] d, final int i, final int j) {
        return 0.5 * d[i][j];
    }

    @Override
    public double dnewfn(final double[][] d, final int k, final int i, final int j, final int ni, final int nj) {
        return (ni * d[i][k] + nj * d[j][k]) / (ni + nj);
    } // Distance is weighted.

    @Override
    public void drootbranchfn(final double[][] d, final int i, final int j, final int ni, final int nj,
            final $double bi, final $double bj) {
        $(bi, $(bj, 0.5 * d[i][j]));
    }

    public Phylo_upgma(final double[][] dist) throws InstantiationException, IllegalAccessException {
        super(dist);
        makethetree(dist);
    } // This call actually makes the

}
