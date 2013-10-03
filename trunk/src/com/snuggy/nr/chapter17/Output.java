
package com.snuggy.nr.chapter17;

import static com.snuggy.nr.util.Static.*;

import com.snuggy.nr.util.*;

public class Output {

    // Structure for output from ODE solver such as odeint.
    private int kmax; // Current capacity of storage arrays.
    private int nvar;
    private int nsave; // Number of intervals to save at for dense output.
    private boolean dense; // true if dense output requested.
    private int count; // Number of values actually saved.
    private double x1, x2, xout, dxout;
    private double[] xsave; // Results stored in the vector xsave[0..count-1]
                            // and the
    private double[][] ysave; // matrix ysave[0..nvar-1][0..count-1].

    public Output() {
        kmax = (-1);
        dense = (false);
        count = (0);
    }

    public int kmax() {
        return kmax;
    }

    public boolean dense() {
        return dense;
    }

    public double[] xsave() {
        return xsave;
    }
    
    public double[][] ysave() {
        return ysave;
    }

    public int count() {
        return count;
    }

    // Default constructor gives no output.

    public Output(final int nsavee) {
        // Constructor provides dense output at nsave equally spaced intervals.
        // If nsave  0, output is saved only at the actual integration steps.
        kmax = (500);
        nsave = (nsavee);
        count = (0);
        xsave = doub_arr(kmax);
        dense = nsave > 0 ? true : false;
    }

    public void init(final int neqn, final double xlo, final double xhi) {
        // Called by Odeint constructor, which passes neqn, the number of
        // equations, xlo, the starting point of the integration, and xhi,
        // the ending point.
        nvar = neqn;
        if (kmax == -1)
            return;
        ysave = doub_mat(nvar, kmax);
        if (dense) {
            x1 = xlo;
            x2 = xhi;
            xout = x1;
            dxout = (x2 - x1) / nsave;
        }
    }

    public void resize() {
        // Resize storage arrays by a factor of two, keeping saved data.
        int kold = kmax;
        kmax *= 2;
        double[] tempvec = doub_arr(xsave);
        xsave = doub_arr(kmax);
        for (int k = 0; k < kold; k++)
            xsave[k] = tempvec[k];
        double[][] tempmat = doub_mat(ysave);
        ysave = doub_mat(nvar, kmax);
        for (int i = 0; i < nvar; i++)
            for (int k = 0; k < kold; k++)
                ysave[i][k] = tempmat[i][k];
    }

    public void save_dense(final IStepperBS s, final double xout, final double h) throws NRException {
        // Invokes dense_out function of stepper routine to produce output at
        // xout. Normally called by out rather than directly. Assumes that xout
        // is between xold and xold+h, where the stepper must keep track of
        // xold, the location of the previous step, and x=xold+h, the current
        // step.
        if (count == kmax)
            resize();
        for (int i = 0; i < nvar; i++)
            ysave[i][count] = s.dense_out(i, xout, h);
        xsave[count++] = xout;
    }

    public void save(final double x, final double[] y) {
        // Saves values of current x and y.
        if (kmax <= 0)
            return;
        if (count == kmax)
            resize();
        for (int i = 0; i < nvar; i++)
            ysave[i][count] = y[i];
        xsave[count++] = x;
    }

    public void out(final int nstp, final double x, final double[] y, final IStepperBS s, final double h)
            throws NRException {
        // Typically called by Odeint to produce dense output. Input variables
        // are nstp, the current step number, the current values of x and y,
        // the stepper s, and the stepsize h. A call with nstp=-1 saves the
        // initial values. The routine checks whether x is greater than the
        // desired output point xout. If so, it calls save_dense.
        if (!dense)
            throw new NRException("dense output not set in Output!");
        if (nstp == -1) {
            save(x, y);
            xout += dxout;
        } else {
            while ((x - xout) * (x2 - x1) > 0.0) {
                save_dense(s, xout, h);
                xout += dxout;
            }
        }
    }

}
