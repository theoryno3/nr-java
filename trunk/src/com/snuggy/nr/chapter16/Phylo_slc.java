
package com.snuggy.nr.chapter16;

import static com.snuggy.nr.refs.Refs.*;
import static com.snuggy.nr.util.Static.*;

import com.snuggy.nr.refs.*;

public class Phylo_slc extends Phylagglom {

    // Derived class implementing the single linkage clustering method.
    @Override
    public void premin(final double[][] d, final int[] nextp) {
    } // No pre-min calculations.

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
        return MIN(d[i][k], d[j][k]);
    } // New-node distance is min of children.

    @Override
    public void drootbranchfn(final double[][] d, final int i, final int j, final int ni, final int nj,
            final $double bi, final $double bj) {
        $(bi, $(bj, 0.5 * d[i][j]));
    }

    public Phylo_slc(final double[][] dist) throws InstantiationException, IllegalAccessException {
        super(dist);
        makethetree(dist);
    } // This call actually makes the tree.

}
