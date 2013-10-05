
package com.snuggy.nr.chapter17;

import static com.snuggy.nr.refs.Refs.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.chapter02.*;
import com.snuggy.nr.refs.*;
import com.snuggy.nr.util.*;

public class StepperSie extends StepperBS {

    // Semi-implicit extrapolation step for integrating sti ODEs, with
    // monitoring of local truncation error to adjust stepsize.
    // typedef D Dtype;
    private static final int KMAXX = 12, IMAXX = KMAXX + 1;
    // KMAXX is the maximum number of rows used in the extrapolation.
    private int k_targ; // Optimal row number for convergence.
    private final int[] nseq; // Stepsize sequence.
    private final double[] cost; // Ak.
    private final double[][] table; // Extrapolation tableau.
    private final double[][] dfdy; // f 0
    private final double[] dfdx; // @f=@x (for compatibility with StepperRoss; not
                           // used.)
    private double jac_redo; // Criterion for recomputing Jacobian.
    private boolean calcjac; // True if Jacobian is current.
    private double theta; // Recompute Jacobian if theta > jac_redo.
    private final double[][] a;
    private int kright; // Used in dense output.
    private final double[][] coeff; // Coecients in extrapolation tableau.
    private final double[][] fsave; // Stores right-hand sides for dense output.
    private final double[] dens; // Stores quantities for dense interpolating
                           // polynomial.
    private final double[] factrl; // Factorials.

    private static final double costfunc = 1.0, costjac = 5.0, costlu = 1.0, costsolve = 1.0;

    // StepperSie(VecDoub_IO &yy, VecDoub_IO &dydxx, Doub &xx, const Doub atol,
    // const Doub rtol, bool dens);
    // void step(const Doub htry,D &derivs);
    // bool dy(VecDoub_I &y, const Doub htot, const Int k, VecDoub_O &yend,
    // Int &ipt,VecDoub_I &scale,D &derivs);
    // void polyextr(const Int k, MatDoub_IO &table, VecDoub_IO &last);
    // void prepare_dense(const Doub h,VecDoub_I &ysav,VecDoub_I &scale,
    // const Int k, Doub &error);
    // Doub dense_out(const Int i,const Doub x,const Doub h);
    // void dense_interp(const Int n, VecDoub_IO &y, const Int imit);

    public StepperSie(final double[] yy, final double[] dydxx, final $double xx, final double atoll, final double rtoll,
            boolean densflag) throws NRException {
        super(yy, dydxx, xx, atoll, rtoll, densflag);
        nseq = int_vec(IMAXX);
        cost = doub_vec(IMAXX);
        table = doub_mat(KMAXX, n);
        dfdy = doub_mat(n, n);
        dfdx = doub_vec(n);
        calcjac = (false);
        a = doub_mat(n, n);
        coeff = doub_mat(IMAXX, IMAXX);
        fsave = doub_mat((IMAXX - 1) * (IMAXX + 1) / 2 + 2, n);
        dens = doub_vec((IMAXX + 2) * n);
        factrl = doub_vec(IMAXX);
        // Input to the constructor are the dependent variable y[0..n-1] and
        // its derivative dydx[0..n-1] at the starting value of the independent
        // variable x. Also input are the absolute and relative tolerances,
        // atol and rtol, and the boolean dense, which is true if dense output
        // is required.
        // The cost of a Jacobian is taken to be 5 function evaluations.
        // Performance is not too sensitive to the value used.
        EPS = EPS(); // numeric_limits<Doub>::epsilon();
        jac_redo = MIN(1.0e-4, rtol);
        theta = 2.0 * jac_redo; // Make sure Jacobian is computed on rst step.
        nseq[0] = 2; // Sequence is dierent from StepperBS.
        nseq[1] = 3;
        for (int i = 2; i < IMAXX; i++)
            nseq[i] = 2 * nseq[i - 2];
        cost[0] = costjac + costlu + nseq[0] * (costfunc + costsolve);
        for (int k = 0; k < KMAXX; k++)
            cost[k + 1] = cost[k] + (nseq[k + 1] - 1) * (costfunc + costsolve) + costlu;
        hnext = -1.0e99; // Impossible value.
        double logfact = -log10(rtol + atol) * 0.6 + 0.5;
        k_targ = MAX(1, MIN(KMAXX - 1, Int(logfact))); // Initial estimate of
                                                       // optimal k.
        for (int k = 0; k < IMAXX; k++) { // Coecients in equation (17.3.8),
                                          // but ra
            for (int l = 0; l < k; l++) { // tio not squared.
                double ratio = Doub(nseq[k]) / nseq[l];
                coeff[k][l] = 1.0 / (ratio - 1.0);
            }
        }
        factrl[0] = 1.0;
        for (int k = 0; k < IMAXX - 1; k++)
            factrl[k + 1] = (k + 1) * factrl[k];
    }

    private final double STEPFAC1 = 0.6, STEPFAC2 = 0.93, STEPFAC3 = 0.1, STEPFAC4 = 4.0, STEPFAC5 = 0.5, KFAC1 = 0.7,
            KFAC2 = 0.9;
    private static boolean first_step = true, last_step = false;
    private static boolean forward, reject = false, prev_reject = false;
    private static double errold;

    public void step(final double htry, final Dtype derivs) throws NRException {
        // Attempts a step with stepsize htry. On output, y and x are replaced
        // by their new values, hdid is the stepsize that was actually
        // accomplished, and hnext is the estimated next stepsize.
        // Stepsize and order control parameters are dierent from StepperBS.
        int i, k = 0;
        double fac, h, hnew;
        $double err = $(0.0);
        boolean firstk;
        final double[] hopt = doub_vec(IMAXX), work = doub_vec(IMAXX);
        final double[] ysav = doub_vec(n), yseq = doub_vec(n);
        @SuppressWarnings("unused")
        final double[] ymid = doub_vec(n);
        final double[] scale = doub_vec(n);
        work[0] = 1.e30;
        h = htry;
        forward = h > 0 ? true : false;
        for (i = 0; i < n; i++)
            ysav[i] = y[i]; // Save the starting values.
        if (h != hnext && !first_step) { // h gets reset in Odeint for the last
                                         // step.
            last_step = true;
        }
        if (reject) { // Previous step was rejected.
            prev_reject = true;
            last_step = false;
            theta = 2.0 * jac_redo; // Make sure Jacobian gets recomputed.
        }
        for (i = 0; i < n; i++)
            // Initial scaling.
            scale[i] = atol + rtol * abs(y[i]);
        reject = false;
        firstk = true;
        hnew = abs(h);
        compute_jac: // Restart here if Jacobian error too big.
        while (true) {
	        if (theta > jac_redo && !calcjac) { // Evaluate Jacobian.
	            derivs.jacobian(x.$(), y, dfdx, dfdy);
	            calcjac = true;
	        }
	        while (firstk || reject) { // Loop until step accepted.
	            h = forward ? hnew : -hnew;
	            firstk = false;
	            reject = false;
	            if (abs(h) <= abs(x.$()) * EPS)
	                throw new NRException("step size underflow in StepperSie");
	            $int ipt = $(-1); // Initialize counter for saving stu.
	            for (k = 0; k <= k_targ + 1; k++) { // The sequence of semi-implicit
	                                                // Euler steps.
	                boolean success = dy(ysav, h, k, yseq, ipt, scale, derivs);
	                if (!success) { // Stability problems, reduce stepsize.
	                    reject = true;
	                    hnew = abs(h) * STEPFAC5;
	                    break;
	                }
	                if (k == 0)
	                    $$(y, yseq);
	                else
	                    // Store result in tableau.
	                    for (i = 0; i < n; i++)
	                        table[k - 1][i] = yseq[i];
	                if (k != 0) {
	                    polyextr(k, table, y); // Perform extrapolation.
	                    $(err, 0.0); // Compute normalized error estimate errk.
	                    for (i = 0; i < n; i++) {
	                        scale[i] = atol + rtol * abs(ysav[i]);
	                        err.$(err.$() + SQR((y[i] - table[0][i]) / scale[i]));
	                    }
	                    err.$(sqrt(err.$() / n));
	                    if (err.$() > 1.0 / EPS || (k > 1 && err.$() >= errold)) {
	                        reject = true; // Stability problems, reduce stepsize.
	                        hnew = abs(h) * STEPFAC5;
	                        break;
	                    }
	                    errold = max(4.0 * err.$(), 1.0);
	                    double expo = 1.0 / (k + 1);
	                    // Compute optimal stepsize for this order. Note k instead
	                    // of 2k in exponent.
	                    double facmin = pow(STEPFAC3, expo);
	                    if (err.$() == 0.0)
	                        fac = 1.0 / facmin;
	                    else {
	                        fac = STEPFAC2 / pow(err.$() / STEPFAC1, expo);
	                        fac = MAX(facmin / STEPFAC4, MIN(1.0 / facmin, fac));
	                    }
	                    hopt[k] = abs(h * fac);
	                    work[k] = cost[k] / hopt[k]; // Work per unit step
	                                                 // (17.3.13).
	                    if ((first_step || last_step) && err.$() <= 1.0)
	                        break;
	                    if (k == k_targ - 1 && !prev_reject && !first_step && !last_step) {
	                        if (err.$() <= 1.0) // Converged within order window.
	                            break;
	                        else if (err.$() > nseq[k_targ] * nseq[k_targ + 1] * 4.0) {
	                            reject = true; // No convergence expected by
	                                           // k_targ+1.
	                            k_targ = k;
	                            if (k_targ > 1 && work[k - 1] < KFAC1 * work[k])
	                                k_targ--;
	                            hnew = hopt[k_targ];
	                            break;
	                        }
	                    }
	                    if (k == k_targ) {
	                        if (err.$() <= 1.0) // Converged within order window.
	                            break;
	                        else if (err.$() > nseq[k + 1] * 2.0) {
	                            reject = true; // No convergence expected by
	                                           // k_targ+1.
	                            if (k_targ > 1 && work[k - 1] < KFAC1 * work[k])
	                                k_targ--;
	                            hnew = hopt[k_targ];
	                            break;
	                        }
	                    }
	                    if (k == k_targ + 1) {
	                        if (err.$() > 1.0) {
	                            reject = true;
	                            if (k_targ > 1 && work[k_targ - 1] < KFAC1 * work[k_targ])
	                                k_targ--;
	                            hnew = hopt[k_targ];
	                        }
	                        break;
	                    }
	                }
	            } // Go back and try next k.
	            if (reject) { // Arrive here from any break in for loop.
	                prev_reject = true;
	                if (!calcjac) {
	                    theta = 2.0 * jac_redo;
	                    continue compute_jac;
	                }
	            }
	        } // Go back if step was rejected.
	        break;
        }
        calcjac = false; // Successful step. Allow Jacobian to be re
        if (dense) // computed if theta too big.
            prepare_dense(h, ysav, scale, k, err);
        xold = x.$(); // Used by dense output.
        x.$(x.$() + h);
        hdid = h;
        first_step = false;
        int kopt; // Determine optimal order for next step.
        if (k == 1)
            kopt = 2;
        else if (k <= k_targ) {
            kopt = k;
            if (work[k - 1] < KFAC1 * work[k])
                kopt = k - 1;
            else if (work[k] < KFAC2 * work[k - 1])
                kopt = MIN(k + 1, KMAXX - 1);
        } else {
            kopt = k - 1;
            if (k > 2 && work[k - 2] < KFAC1 * work[k - 1])
                kopt = k - 2;
            if (work[k] < KFAC2 * work[kopt])
                kopt = MIN(k, KMAXX - 1);
        }
        if (prev_reject) { // After a rejected step neither order nor step
            k_targ = MIN(kopt, k); // size should increase.
            hnew = MIN(abs(h), hopt[k_targ]);
            prev_reject = false;
        } else { // Stepsize control for next step.
            if (kopt <= k)
                hnew = hopt[kopt];
            else {
                if (k < k_targ && work[k] < KFAC2 * work[k - 1])
                    hnew = hopt[k] * cost[kopt + 1] / cost[k];
                else
                    hnew = hopt[k] * cost[kopt] / cost[k];
            }
            k_targ = kopt;
        }
        if (forward)
            hnext = hnew;
        else
            hnext = -hnew;
    }

    public boolean dy(final double[] y, final double htot, final int k, final double[] yend, final $int ipt,
            final double[] scale, final Dtype derivs) throws NRException {
        // Semi-implicit Euler step. Inputs are y, H, k and scale[0..n-1]. The
        // output is returned as yend[0..n-1]. The counter ipt keeps track of
        // saving the right-hand sides in the correct locations for dense
        // output.
        final double[] del = doub_vec(n), ytemp = doub_vec(n), dytemp = doub_vec(n);
        int nstep = nseq[k];
        double h = htot / nstep; // Stepsize this trip.
        for (int i = 0; i < n; i++) { // Set up the matrix 1=h f 0.
            for (int j = 0; j < n; j++)
                a[i][j] = -dfdy[i][j];
            a[i][i] += 1.0 / h;
        }
        LUdcmp alu = new LUdcmp(a); // LU decomposition of the matrix.
        double xnew = x.$() + h; // Special step for nonautonomous system.
        derivs.eval(xnew, y, del);
        for (int i = 0; i < n; i++)
            ytemp[i] = y[i];
        alu.solve(del, del);
        if (dense && nstep == k + 1) {
            ipt.$(ipt.$() + 1);
            for (int i = 0; i < n; i++)
                fsave[ipt.$()][i] = del[i];
        }
        for (int nn = 1; nn < nstep; nn++) { // General step.
            for (int i = 0; i < n; i++)
                ytemp[i] += del[i];
            xnew += h;
            derivs.eval(xnew, ytemp, yend);
            if (nn == 1 && k <= 1) { // Stability test and test for recomputing
                                     // Jaco-
                double del1 = 0.0; // bian.
                for (int i = 0; i < n; i++)
                    del1 += SQR(del[i] / scale[i]);
                del1 = sqrt(del1);
                derivs.eval(x.$() + h, ytemp, dytemp);
                for (int i = 0; i < n; i++)
                    del[i] = dytemp[i] - del[i] / h;
                alu.solve(del, del);
                double del2 = 0.0;
                for (int i = 0; i < n; i++)
                    del2 += SQR(del[i] / scale[i]);
                del2 = sqrt(del2);
                theta = del2 / MAX(1.0, del1);
                if (theta > 1.0)
                    return false;
            }
            alu.solve(yend, del);
            if (dense && nn >= nstep - k - 1) {
                ipt.$(ipt.$() + 1);
                for (int i = 0; i < n; i++)
                    fsave[ipt.$()][i] = del[i];
            }
        }
        for (int i = 0; i < n; i++)
            // Last step.
            yend[i] = ytemp[i] + del[i];
        return true;
    }

    public void polyextr(final int k, final double[][] table, final double[] last) {
        // Use polynomial extrapolation to evaluate l functions at h D 0. This
        // routine is identical to the routine in StepperBS.
        int l = last.length;
        for (int j = k - 1; j > 0; j--)
            for (int i = 0; i < l; i++)
                table[j - 1][i] = table[j][i] + coeff[k][j] * (table[j][i] - table[j - 1][i]);
        for (int i = 0; i < l; i++)
            last[i] = table[0][i] + coeff[k][0] * (table[0][i] - last[i]);
    }

    public void prepare_dense(final double h, final double[] ysav, final double[] scale, final int k,
            final $double error) {
        // Store coecients of interpolating polynomial for dense output in dens
        // array. Input stepsize h, function at beginning of interval
        // ysav[0..n-1], scale factor atolCjyjrtol in scale[0..n-1], and column
        // k in which convergence was achieved. Output interpolation error in
        // error.
        kright = k;
        for (int i = 0; i < n; i++) {
            dens[i] = ysav[i];
            dens[n + i] = y[i];
        }
        for (int klr = 0; klr < kright; klr++) { // Compute dierences.
            if (klr >= 1) {
                for (int kk = klr; kk <= k; kk++) {
                    int lbeg = ((kk + 3) * kk) / 2;
                    int lend = lbeg - kk + 1;
                    for (int l = lbeg; l >= lend; l--)
                        for (int i = 0; i < n; i++)
                            fsave[l][i] = fsave[l][i] - fsave[l - 1][i];
                }
            }
            for (int kk = klr; kk <= k; kk++) { // Compute derivatives at right
                                                // end.
                double facnj = nseq[kk];
                facnj = pow(facnj, klr + 1) / factrl[klr + 1];
                int ipt = ((kk + 3) * kk) / 2;
                int krn = (kk + 2) * n;
                for (int i = 0; i < n; i++) {
                    dens[krn + i] = fsave[ipt][i] * facnj;
                }
            }
            for (int j = klr + 1; j <= k; j++) {
                double dblenj = nseq[j];
                for (int l = j; l >= klr + 1; l--) {
                    double factor = dblenj / nseq[l - 1] - 1.0;
                    for (int i = 0; i < n; i++) {
                        int krn = (l + 2) * n + i;
                        dens[krn - n] = dens[krn] + (dens[krn] - dens[krn - n]) / factor;
                    }
                }
            }
        }
        for (int in = 0; in < n; in++) { // Compute coecients of the
                                         // interpolation poly
            for (int j = 1; j <= kright + 1; j++) { // nomial.
                int ii = n * j + in;
                dens[ii] = dens[ii] - dens[ii - n];
            }
        }
    }

    public double dense_out(final int i, final double x, final double h) {
        // Evaluate interpolating polynomial for y[i] at location x, where xold
        //  x  xoldCh.
        double theta = (x - xold) / h;
        int k = kright;
        double yinterp = dens[(k + 1) * n + i];
        for (int j = 1; j <= k; j++)
            yinterp = dens[(k + 1 - j) * n + i] + yinterp * (theta - 1.0);
        return dens[i] + yinterp * theta;
    }
}
