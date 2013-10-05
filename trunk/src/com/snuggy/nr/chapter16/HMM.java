
package com.snuggy.nr.chapter16;

import static com.snuggy.nr.refs.Refs.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.refs.*;
import com.snuggy.nr.util.*;

@Deprecated @Broken
public class HMM {

    // Structure for a hidden Markov model and its methods.
    private final $$double2d a, b; // Transition matrix and symbol probability matrix.
    private final $$int1d obs; // Observed data.
    private int fbdone;
    private int mstat, nobs, ksym; // Number of states, observations, and
                                   // symbols.
    private int lrnrm;
    private final $$double2d alpha, beta, pstate; // Matrices ?, ?, and Pi .t/.
    private final $$int1d arnrm, brnrm;
    private double BIG, BIGI, lhood;

    // HMM(MatDoub_I &aa, MatDoub_I &bb, VecInt_I &obs); Constructor; see below.
    // void forwardbackward(); HMM state estimation.
    // void baumwelch(); HMM parameter re-estimation.

    public double loglikelihood() {
        // Returns the log-likelihood computed by forwardbackward().
        return log(lhood) + lrnrm * log(BIGI);
    }
    
    public double[][] pstate() {
        return pstate.$();
    }
    
    public HMM(final double[][] aa, final double[][] bb, final int[] obss) throws NRException {
        // Constructor. Input are the transition matrix aa, the symbol
        // probability matrix bb, and the observed vector of symbols obss.
        // Local copies are made, so the input quantities need not be preserved
        // by the calling program.
        a = $$(aa);
        b = $$(bb);
        obs = $$(obss);
        fbdone = (0);
        mstat = (nrows(a.$()));
        nobs = (obs.$().length);
        ksym = (ncols(b.$()));
        alpha = $$(doub_mat(nobs, mstat));
        beta = $$(doub_mat(nobs, mstat));
        pstate = $$(doub_mat(nobs, mstat));
        arnrm = $$(int_vec(nobs));
        brnrm = $$(int_vec(nobs));
        BIG = (1.e20);
        BIGI = (1. / BIG);
        int i, j, k;
        double sum;
        // Although space constraints make us generally stingy about printing
        // code for checking input, we will save you a lot of grief by doing
        // so in this case. If you get “matrix not normalized” errors, you
        // probably have your matrix transposed. Note that normalization
        // errors <1% are silently fixed.
        if (ncols(a.$()) != mstat)
            throw new NRException("transition matrix not square");
        if (nrows(b.$()) != mstat)
            throw new NRException("symbol prob matrix wrong size");
        for (i = 0; i < nobs; i++) {
            if (obs.$()[i] < 0 || obs.$()[i] >= ksym)
                throw new NRException("bad data in obs");
        }
        for (i = 0; i < mstat; i++) {
            sum = 0.;
            for (j = 0; j < mstat; j++)
                sum += a.$()[i][j];
            if (abs(sum - 1.) > 0.01)
                throw new NRException("transition matrix not normalized");
            for (j = 0; j < mstat; j++)
                a.$()[i][j] /= sum;
        }
        for (i = 0; i < mstat; i++) {
            sum = 0.;
            for (k = 0; k < ksym; k++)
                sum += b.$()[i][k];
            if (abs(sum - 1.) > 0.01)
                throw new NRException("symbol prob matrix not normalized");
            for (k = 0; k < ksym; k++)
                b.$()[i][k] /= sum;
        }
    }

    // Now, to actually do the forward-backward estimation, you call the
    // function forwardbackward. This fills the matrix pstate, so that
    // pstateti D Pt .i /. It also sets the internal variables lhood and
    // lrnrm so that the function loglikelihood returns the logarithm of L.
    // Don’t be surprised at how large in magnitude this (negative) number
    // can be. The probability of any particular data set of more than
    // trivial length is astronomically small! In the following code, the
    // quantities BIG, BIGI, arnrm, brnrm, and lrnrm all relate to dealing
    // with values that would far underflow an ordinary floating format. The
    // basic idea is to renormalize as necessary, keeping track of the
    // accumulated number of renormalizations. At the end, when an ?, a ?,
    // and an L are combined, probability values of reasonable magnitude result.

    public void forwardbackward() {
        // HMM forward-backward algorithm. Using the stored a, b, and obs
        // matrices, the matrices alpha, beta, and pstate are calculated.
        // The latter is the state estimation of the model, given the data.
        int i, j, t;
        double sum, asum, bsum;
        for (i = 0; i < mstat; i++)
            alpha.$()[0][i] = b.$()[i][obs.$()[0]];
        arnrm.$()[0] = 0;
        for (t = 1; t < nobs; t++) { // Forward pass.
            asum = 0;
            for (j = 0; j < mstat; j++) {
                sum = 0.;
                for (i = 0; i < mstat; i++)
                    sum += alpha.$()[t - 1][i] * a.$()[i][j] * b.$()[j][obs.$()[t]];
                alpha.$()[t][j] = sum;
                asum += sum;
            }
            arnrm.$()[t] = arnrm.$()[t - 1]; // Renormalize the ?’s as necessary to
                                     // avoid
            // underflow, keeping track of how many renormalizations for each ?.
            if (asum < BIGI) {
                ++arnrm.$()[t];
                for (j = 0; j < mstat; j++)
                    alpha.$()[t][j] *= BIG;
            }
        }
        for (i = 0; i < mstat; i++)
            beta.$()[nobs - 1][i] = 1.;
        brnrm.$()[nobs - 1] = 0;
        for (t = nobs - 2; t >= 0; t--) { // Backward pass.
            bsum = 0.;
            for (i = 0; i < mstat; i++) {
                sum = 0.;
                for (j = 0; j < mstat; j++)
                    sum += a.$()[i][j] * b.$()[j][obs.$()[t + 1]] * beta.$()[t + 1][j];
                beta.$()[t][i] = sum;
                bsum += sum;
            }
            brnrm.$()[t] = brnrm.$()[t + 1];
            if (bsum < BIGI) { // Similarly, renormalize the ?’s as necessary.
                ++brnrm.$()[t];
                for (j = 0; j < mstat; j++)
                    beta.$()[t][j] *= BIG;
            }
        }
        lhood = 0.; // Overall likelihood is lhood with lnorm renormal
        for (i = 0; i < mstat; i++)
            lhood += alpha.$()[0][i] * beta.$()[0][i]; // izations.
        lrnrm = arnrm.$()[0] + brnrm.$()[0];
        while (lhood < BIGI) {
            lhood *= BIG;
            lrnrm++;
        }
        for (t = 0; t < nobs; t++) { // Get state probabilities from ?’s and
                                     // ?’s.
            sum = 0.;
            for (i = 0; i < mstat; i++)
                sum += (pstate.$()[t][i] = alpha.$()[t][i] * beta.$()[t][i]);
            // The next line is an equivalent calculation of sum. But we’d
            // rather
            // have the normalization of the Pi .t/’s be more immune to roundoff
            // error. Hence we do the above sum for each value of t.
            // // sum = lhood*pow(BIGI, lrnrm - arnrm[t] - brnrm[t]);
            for (i = 0; i < mstat; i++)
                pstate.$()[t][i] /= sum;
        }
        fbdone = 1; // Flag prevents misuse of baumwelch(), later.
    }

    public void baumwelch() throws NRException {
        // Baum-Welch re-estimation of the stored matrices a and b, using the
        // data obs and the matrices alpha and beta as computed by
        // forwardbackward() (which must be called first). The previous values
        // of a and b are overwritten.
        int i, j, k, t;
        double num, denom, term;
        final double[][] bnew = doub_mat(mstat, ksym);
        final int[] powtab = int_vec(10); // Fill table of powers of BIGI.
        for (i = 0; i < 10; i++)
            powtab[i] = (int) pow(BIGI, i - 6);
        if (fbdone != 1)
            throw new NRException("must do forwardbackward first");
        for (i = 0; i < mstat; i++) { // Looping over i, get denominators and
                                      // new b.
            denom = 0.;
            for (k = 0; k < ksym; k++)
                bnew[i][k] = 0.;
            for (t = 0; t < nobs - 1; t++) {
                term = (alpha.$()[t][i] * beta.$()[t][i] / lhood) * powtab[arnrm.$()[t] + brnrm.$()[t] - lrnrm + 6];
                denom += term;
                bnew[i][obs.$()[t]] += term;
            }
            for (j = 0; j < mstat; j++) { // Inner loop over j gets elements of
                                          // a.
                num = 0.;
                for (t = 0; t < nobs - 1; t++) {
                    num += alpha.$()[t][i] * b.$()[j][obs.$()[t + 1]] * beta.$()[t + 1][j]
                            * powtab[arnrm.$()[t] + brnrm.$()[t + 1] - lrnrm + 6] / lhood;
                }
                a.$()[i][j] *= (num / denom);
            }
            for (k = 0; k < ksym; k++)
                bnew[i][k] /= denom;
        }
        $$(b, bnew);
        fbdone = 0; // Don’t let this routine be called again until forward}
        // backward() has been called.

    }

}
