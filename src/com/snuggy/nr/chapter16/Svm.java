
package com.snuggy.nr.chapter16;

import static com.snuggy.nr.refs.Refs.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.chapter07.*;
import com.snuggy.nr.chapter08.*;
import com.snuggy.nr.refs.*;
import com.snuggy.nr.util.*;

public class Svm {

    // Class for solving SVM problems by the SOR method.
    private Svmgenkernel gker; // Reference bound to user’s kernel (and data).
    private int m, fnz, fub, niter;
    private $$double1d alph, alphold; // Vectors of ?’s before and after a step.
    private Ran ran; // Random number generator.
    private boolean alphinit;
    private double dalph; // Change in norm of the ?’s in one step.

    public Svm(final Svmgenkernel inker) throws NRException {
        gker = (inker);
        m = (gker.y().length);
        alph = $$(doub_vec(m));
        alphold = $$(new double[m]);
        ran = new Ran(21);
        alphinit = (false);
    }

    // Constructor binds the user’s kernel and allocates storage.

    public double relax(final double lambda, final double om) throws NRException {
        // Perform one group of relaxation steps: a single step over all
        // the ?’s, and multiple steps over only the interior ?’s.
        int iter, j, jj, k, kk;
        double sum; // Index when ?’s are sorted by value.
        final double[] pinsum = doub_vec(m); // Stored sums over noninterior
                                        // variables.
        if (alphinit == false) { // Start all ?’s at 0.
            for (j = 0; j < m; j++)
                alph.$()[j] = 0.;
            alphinit = true;
        }
        $$(alphold, alph); // Save old ?’s.
        // Here begins the relaxation pass over all the ?’s.
        DoubleIndexx x = new DoubleIndexx(alph.$()); // Sort ?’s, then find first
                                                 // nonzero one.
        for (fnz = 0; fnz < m; fnz++)
            if (alph.$()[x.indx()[fnz]] != 0.)
                break;
        for (j = fnz; j < m - 2; j++) { // Randomly permute all the nonzero ?’s.
            //k = j + (ran.int32() % (m - j));
            k = j + umod(ran.int32(), (m - j));
            SWAP(x.indx(), j, k);
        }
        for (jj = 0; jj < m; jj++) { // Main loop over ?’s.
            j = x.indx()[jj];
            sum = 0.;
            for (kk = fnz; kk < m; kk++) { // Sums start with first nonzero.
                k = x.indx()[kk];
                sum += (gker.ker()[j][k] + 1.) * gker.y()[k] * alph.$()[k];
            }
            alph.$()[j] = alph.$()[j] - (om / (gker.ker()[j][j] + 1.)) * (gker.y()[j] * sum - 1.);
            alph.$()[j] = MAX(0., MIN(lambda, alph.$()[j])); // Projection operator.
            if ((jj < fnz) && (alph.$()[j] != 0))
                SWAP(x.indx(), --fnz, jj);
        } // (Above) Make an ? active if it becomes nonzero.
          // Here begins the relaxation passes over the interior ?’s.
        DoubleIndexx y = new DoubleIndexx(alph.$()); // Sort. Identify interior ?’s.
        for (fnz = 0; fnz < m; fnz++)
            if (alph.$()[y.indx()[fnz]] != 0.)
                break;
        for (fub = fnz; fub < m; fub++)
            if (alph.$()[y.indx()[fub]] == lambda)
                break;
        for (j = fnz; j < fub - 2; j++) { // Permute.
            k = j + umod(ran.int32(), (fub - j));
            SWAP(y.indx(), j, k);
        }
        for (jj = fnz; jj < fub; jj++) { // Compute sums over pinned ?’s just
            j = y.indx()[jj]; // once.
            sum = 0.;
            for (kk = fub; kk < m; kk++) {
                k = y.indx()[kk];
                sum += (gker.ker()[j][k] + 1.) * gker.y()[k] * alph.$()[k];
            }
            pinsum[jj] = sum;
        }
        niter = MAX(Int(0.5 * (m + 1.0) * (m - fnz + 1.0) / (SQR(fub - fnz + 1.0))), 1);
        // Calculate a number of iterations that will take about half as
        // long as the full pass just completed.
        for (iter = 0; iter < niter; iter++) { // Main loop over ?’s.
            for (jj = fnz; jj < fub; jj++) {
                j = y.indx()[jj];
                sum = pinsum[jj];
                for (kk = fnz; kk < fub; kk++) {
                    k = y.indx()[kk];
                    sum += (gker.ker()[j][k] + 1.) * gker.y()[k] * alph.$()[k];
                }
                alph.$()[j] = alph.$()[j] - (om / (gker.ker()[j][j] + 1.)) * (gker.y()[j] * sum - 1.);
                alph.$()[j] = MAX(0., MIN(lambda, alph.$()[j]));
            }
        }
        dalph = 0.; // Return change in norm of ? vector.
        for (j = 0; j < m; j++)
            dalph += SQR(alph.$()[j] - alphold.$()[j]);
        return sqrt(dalph);
    }

    public double predict(final int k) {
        // Call only after convergence via repeated calls to relax.
        // Returns the decision rule f.x/ for data point k.
        double sum = 0.;
        for (int j = 0; j < m; j++)
            sum += alph.$()[j] * gker.y()[j] * (gker.ker()[j][k] + 1.0);
        return sum;
    }

    public double predict(final double x_arr[], final int x_off) {
        // Call only after convergence via repeated calls to relax. Returns
        // the decision rule f.x/ for an arbitrary feature vector.
        double sum = 0.;
        for (int j = 0; j < m; j++)
            sum += alph.$()[j] * gker.y()[j] * (gker.kernel(j, x_arr, x_off) + 1.0);
        return sum;
    }

}
