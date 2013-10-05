
package com.snuggy.nr.chapter02;

import static com.snuggy.nr.refs.Refs.*;

import com.snuggy.nr.refs.*;
import com.snuggy.nr.util.*;

public class NRsparseLinbcg extends Linbcg {

    // Here is an example of a derived class that solvesA x D b for a
    // matrix A in NRsparseMat’s compressed column sparse format. A naive
    // diagonal preconditioner is used.

    private NRsparseMat mat;
    private int n;

    public NRsparseLinbcg(final NRsparseMat matrix) {
        mat = new NRsparseMat(matrix);
        n = (mat.nrows());
    }

    // The constructor just binds a reference to your sparse matrix, making
    // it available to asolve and atimes. To solve for a right-hand side,
    // you call this object’s solve method, as defined in the base class.

    public void atimes(final double[] x, final $$double1d r, final int itrnsp) throws NRException {
        if (itrnsp != 0)
            $$(r, mat.atx(x));
        else
            $$(r, mat.ax(x));
    }

    public void asolve(final double[] b, final double[] x, final int itrnsp) {
        int i, j;
        double diag;
        for (i = 0; i < n; i++) {
            diag = 0.0;
            for (j = mat.col_ptr()[i]; j < mat.col_ptr()[i + 1]; j++)
                if (mat.row_ind()[j] == i) {
                    diag = mat.val()[j];
                    break;
                }
            x[i] = (diag != 0.0 ? b[i] / diag : b[i]);
            // The matrix zA is the diagonal part of A. Since the transpose
            // matrix has the same diagonal, the flag itrnsp is not used in 
            // this example.
        }
    }

}
