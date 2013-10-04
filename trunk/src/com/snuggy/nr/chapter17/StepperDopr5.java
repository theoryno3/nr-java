
package com.snuggy.nr.chapter17;

import static com.snuggy.nr.refs.Refs.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.refs.*;
import com.snuggy.nr.util.*;

@Deprecated @Broken
public class StepperDopr5 extends StepperBS {

    // Dormand-Prince fifth-order Runge-Kutta step with monitoring of local
    // truncation error to ensure accuracy and adjust stepsize.
    // typedef D Dtype; Make the type of derivs available to odeint.
    private double[] k2, k3, k4, k5, k6;
    private double[] rcont1, rcont2, rcont3, rcont4, rcont5;
    private double[] dydxnew;

    // StepperDopr5(VecDoub_IO &yy, VecDoub_IO &dydxx, Doub &xx,
    // const Doub atoll, const Doub rtoll, bool dens);
    // void step(const Doub htry,D &derivs);
    // void dy(const Doub h,D &derivs);
    // void prepare_dense(const Doub h,D &derivs);
    // Doub dense_out(const Int i, const Doub x, const Doub h);
    // Doub error();

    class Controller {
        public double hnext, errold;
        public boolean reject;

        // Controller();
        // bool success(const Doub err, Doub &h);
        // Finally, the controller tests whether err 1 and adjusts the stepsize.
        // The default setting is beta D 0 (no PI control). Set beta to 0.04 or
        // 0.08
        // to turn on PI control.

        public Controller() {
            reject = (false);
            errold = (1.0e-4);
        }

        private static final double beta = 0.0, alpha = 0.2 - beta * 0.75, safe = 0.9, minscale = 0.2, maxscale = 10.0;

        // Step size controller for fifth-order Dormand-Prince method.
        public boolean success(final double err, final $double h) {
            // Returns true if err  1, false otherwise. If step was successful,
            // sets hnext to the estimated optimal stepsize for the next step.
            // If the step failed, reduces h appropriately for another try.
            // Set beta to a nonzero value for PI control. beta D 0:04–0.08 is a
            // good default.
            double scale;
            if (err <= 1.0) { // Step succeeded. Compute hnext.
                if (err == 0.0)
                    scale = maxscale;
                else { // PI control if beta ¤ 0.
                    scale = safe * pow(err, -alpha) * pow(errold, beta);
                    if (scale < minscale)
                        scale = minscale; // Ensure minscale  hnext=h 
                                          // maxscale.
                    if (scale > maxscale)
                        scale = maxscale;
                }
                if (reject) // Don’t let step increase if last one was re
                    hnext = h.$() * MIN(scale, 1.0); // jected.
                else
                    hnext = h.$() * scale;
                errold = MAX(err, 1.0e-4); // Bookkeeping for next call.
                reject = false;
                return true;
            } else { // Truncation error too large, reduce stepsize.
                scale = MAX(safe * pow(err, -alpha), minscale);
                h.$(h.$() * scale);
                reject = true;
                return false;
            }
        }

    }

    Controller con;

    // The constructor simply invokes the base class instructor and initializes
    // variables:

    public StepperDopr5(final double[] yy, final double[] dydxx, final $double xx, final double atoll,
            final double rtoll, final boolean dens) throws NRException {
        super(yy, dydxx, xx, atoll, rtoll, dens);
        k2 = doub_arr(n);
        k3 = doub_arr(n);
        k4 = doub_arr(n);
        k5 = doub_arr(n);
        k6 = doub_arr(n);
        rcont1 = doub_arr(n);
        rcont2 = doub_arr(n);
        rcont3 = doub_arr(n);
        rcont4 = doub_arr(n);
        rcont5 = doub_arr(n);
        dydxnew = doub_arr(n);
        // Input to the constructor are the dependent variable y[0..n-1] and
        // its derivative dydx[0..n-1] at the starting value of the independent
        // variable x. Also input are the absolute and relative tolerances,
        // atol and rtol, and the boolean dense, which is true if dense output
        // is required.
        EPS = EPS(); // numeric_limits<Doub>::epsilon();
    }

    // The step method is the actual stepper. It attempts a step, invokes the
    // controller to decide whether to accept the step or try again with a
    // smaller
    // stepsize, and sets up the coefficients in case dense output is needed
    // between x and x C h.

    public void step(final double htry, final Dtype derivs) throws NRException {
        // Attempts a step with stepsize htry. On output, y and x are replaced
        // by their new values, hdid is the stepsize that was actually
        // accomplished, and hnext is the estimated next stepsize.
        $double h = $(htry); // Set stepsize to the initial trial value.
        for (;;) {
            dy(h.$(), derivs); // Take a step.
            double err = error(); // Evaluate accuracy.
            if (con.success(err, h))
                break; // Step rejected. Try again with reduced h set
            if (abs(h.$()) <= abs(x.$()) * EPS) // by controller.
                throw new NRException("stepsize underflow in StepperDopr5");
        }
        if (dense) // Step succeeded. Compute coefficients for dense
            prepare_dense(h.$(), derivs); // output.
        $$(dydx, dydxnew); // Reuse last derivative evaluation for next step.
        y.$(yout);
        xold = x.$(); // Used for dense output.
        x.$(x.$() + (hdid = h.$()));
        hnext = con.hnext;
    }

    // The algorithm routine dy does the six steps plus the seventh FSAL step,
    // and computes ynC1 and the error .

    private static final double c2 = 0.2, c3 = 0.3, c4 = 0.8, c5 = 8.0 / 9.0, a21 = 0.2, a31 = 3.0 / 40.0,
            a32 = 9.0 / 40.0, a41 = 44.0 / 45.0, a42 = -56.0 / 15.0, a43 = 32.0 / 9.0, a51 = 19372.0 / 6561.0,
            a52 = -25360.0 / 2187.0, a53 = 64448.0 / 6561.0, a54 = -212.0 / 729.0, a61 = 9017.0 / 3168.0,
            a62 = -355.0 / 33.0, a63 = 46732.0 / 5247.0, a64 = 49.0 / 176.0, a65 = -5103.0 / 18656.0,
            a71 = 35.0 / 384.0, a73 = 500.0 / 1113.0, a74 = 125.0 / 192.0, a75 = -2187.0 / 6784.0, a76 = 11.0 / 84.0,
            e1 = 71.0 / 57600.0, e3 = -71.0 / 16695.0, e4 = 71.0 / 1920.0, e5 = -17253.0 / 339200.0, e6 = 22.0 / 525.0,
            e7 = -1.0 / 40.0;

    public void dy(final double h, final Dtype derivs) throws NRException {
        // Given values for n variables y[0..n-1] and their derivatives
        // dydx[0..n-1] known at x, use the fifth-order Dormand-Prince
        // Runge-Kutta method to advance the solution over an interval h and
        // store the incremented variables in yout[0..n-1]. Also store an
        // estimate of the local truncation error in yerr using the
        // embedded fourth-order method.
        double[] ytemp = doub_arr(n);
        int i;
        for (i = 0; i < n; i++)
            // First step.
            ytemp[i] = y.$()[i] + h * a21 * dydx.$()[i];
        derivs.eval(x.$() + c2 * h, ytemp, k2); // Second step.
        for (i = 0; i < n; i++)
            ytemp[i] = y.$()[i] + h * (a31 * dydx.$()[i] + a32 * k2[i]);
        derivs.eval(x.$() + c3 * h, ytemp, k3); // Third step.
        for (i = 0; i < n; i++)
            ytemp[i] = y.$()[i] + h * (a41 * dydx.$()[i] + a42 * k2[i] + a43 * k3[i]);
        derivs.eval(x.$() + c4 * h, ytemp, k4); // Fourth step.
        for (i = 0; i < n; i++)
            ytemp[i] = y.$()[i] + h * (a51 * dydx.$()[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
        derivs.eval(x.$() + c5 * h, ytemp, k5); // Fifth step.
        for (i = 0; i < n; i++)
            ytemp[i] = y.$()[i] + h * (a61 * dydx.$()[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
        double xph = x.$() + h;
        derivs.eval(xph, ytemp, k6); // Sixth step.
        for (i = 0; i < n; i++)
            // Accumulate increments with proper weights.
            yout[i] = y.$()[i] + h * (a71 * dydx.$()[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
        derivs.eval(xph, yout, dydxnew); // Will also be first evaluation for
                                         // next step.
        for (i = 0; i < n; i++) {
            // Estimate error as difference between fourth- and fifth-order
            // methods.
            yerr[i] = h * (e1 * dydx.$()[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * dydxnew[i]);
        }
    }

    // The routine prepare_dense uses the coefficients of [4] to set up the
    // dense output quantities. Our coding of the dense output is closely based
    // on that of the Fortran code DOPRI5 of [5].

    private static final double d1 = -12715105075.0 / 11282082432.0, d3 = 87487479700.0 / 32700410799.0,
            d4 = -10690763975.0 / 1880347072.0, d5 = 701980252875.0 / 199316789632.0, d6 = -1453857185.0 / 822651844.0,
            d7 = 69997945.0 / 29380423.0;

    public void prepare_dense(final double h, final Dtype derivs) {
        // Store coefficients of interpolating polynomial for dense output
        // in rcont1...rcont5.
        @SuppressWarnings("unused")
        double[] ytemp = doub_arr(n);
        for (int i = 0; i < n; i++) {
            rcont1[i] = y.$()[i];
            double ydiff = yout[i] - y.$()[i];
            rcont2[i] = ydiff;
            double bspl = h * dydx.$()[i] - ydiff;
            rcont3[i] = bspl;
            rcont4[i] = ydiff - h * dydxnew[i] - bspl;
            rcont5[i] = h * (d1 * dydx.$()[i] + d3 * k3[i] + d4 * k4[i] + d5 * k5[i] + d6 * k6[i] + d7 * dydxnew[i]);
        }
    }

    // The next routine, dense_out, uses the coefficients stored by the previous
    // routine to evaluate the solution at an arbitrary point.

    public double dense_out(final int i, final double x, final double h) {
        // Evaluate interpolating polynomial for y[i] at location x, where
        // xold  x  xold Ch.
        double s = (x - xold) / h;
        double s1 = 1.0 - s;
        return rcont1[i] + s * (rcont2[i] + s1 * (rcont3[i] + s * (rcont4[i] + s1 * rcont5[i])));
    }

    // The error routine converts  into the scaled quantity err.
    public double error() {
        // Use yerr to compute norm of scaled error estimate. A value less than
        // one means the step was successful.
        double err = 0.0, sk;
        for (int i = 0; i < n; i++) {
            sk = atol + rtol * MAX(abs(y.$()[i]), abs(yout[i]));
            err += SQR(yerr[i] / sk);
        }
        return sqrt(err / n);
    }

    @Override
    public int neqn() {
        return neqn;
    }

    @Override
    public double hdid() {
        return hdid;
    }

    @Override
    public double hnext() {
        return hdid;
    }

}
