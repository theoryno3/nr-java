
package com.snuggy.nr.chapter17;

import static java.lang.Math.*;
import static com.snuggy.nr.util.Static.*;
import static com.snuggy.nr.refs.Refs.*;

import com.snuggy.nr.chapter02.*;
import com.snuggy.nr.refs.*;
import com.snuggy.nr.util.*;

@Deprecated @Broken
public class StepperRoss extends StepperBS {

    // Fourth-order stiffly stable Rosenbrock step for integrating stiff ODEs,
    // with monitoring of local truncation error to adjust stepsize.
    // typedef D Dtype; Make the type of derivs available to odeint.
    private double[][] dfdy; // f 0
    private double[] dfdx; // @f=@x
    private double[] k1, k2, k3, k4, k5, k6;
    private double[] cont1, cont2, cont3, cont4;
    private double[][] a;

    // StepperRoss(double[]_IO &yy, double[]_IO &dydxx, double &xx, final double
    // atoll,
    // final double rtoll, boolean dens);
    // void step(final double htry,D &derivs);
    // void dy(final double h,D &derivs);
    // void prepare_dense(final double h,double[]_I &dydxnew);
    // double dense_out(final int i, final double x, final double h);
    // double error();

    class Controller {
        public double hnext;
        public boolean reject;
        public boolean first_step; // first_step, errold, and hold are used by
        public double errold; // the predictive controller.
        public double hold;

        // Controller();
        // boolean success(double err, double &h);
        public Controller() {
            reject = (false);
            first_step = (true);
        }

        // Step size controller for fourth-order Rosenbrock method.

        private static final double safe = 0.9, fac1 = 5.0, fac2 = 1.0 / 6.0;

        public boolean success(final double err, final $double h) {
            // Returns true if err  1, false otherwise. If step was successful,
            // sets hnext to the estimated optimal stepsize for the next step.
            // If the step failed, reduces h appropriately for another try.
            double fac = MAX(fac2, MIN(fac1, pow(err, 0.25) / safe));
            double hnew = h.$() / fac; // Ensure 1=fac1  hnew=h  1=fac2.
            if (err <= 1.0) { // Step succeeded.
                if (!first_step) { // Predictive control.
                    double facpred = (hold / h.$()) * pow(err * err / errold, 0.25) / safe;
                    facpred = MAX(fac2, MIN(fac1, facpred));
                    fac = MAX(fac, facpred);
                    hnew = h.$() / fac;
                }
                first_step = false;
                hold = h.$();
                errold = MAX(0.01, err);
                if (reject) // Don’t let step increase if last one was rejected.
                    hnew = (h.$() >= 0.0 ? MIN(hnew, h.$()) : MAX(hnew, h.$()));
                hnext = hnew;
                reject = false;
                return true;
            } else { // Truncation error too large, reduce stepsize.
                h.$(hnew);
                reject = true;
                return false;
            }
        }
    }

    Controller con;

    // The implementation will seem very familiar if you’ve looked at
    // StepperDopr5, the explicit Runge-Kutta routine. Note that in the
    // algorithm
    // routine dy of StepperRoss, the linear equations (17.5.32) are solved by
    // first computing the LU decomposition of the matrix 1=h  f 0 using the
    // routine LUdcmp. Then the six gi are found by backsubstitution of the six
    // different right-hand sides using the routine solve in LUdcmp. Thus each
    // step of the integration requires one call to jacobian and six calls to
    // derivs (one call outside dy and five calls inside). The evaluation of
    // the Jacobian matrix is roughly equivalent to N evaluations of the
    // right-hand side f (although it can often be less than this, especially
    // if commonality of code can be exploited). Thus this scheme involves
    // about N C6 function evaluations per step. Note that if N is large and
    // the Jacobian matrix is sparse, you should replace the LU decomposition
    // by a suitable sparse matrix procedure.

    public StepperRoss(final double[] yy, final double[] dydxx, final $double xx, final double atoll,
            final double rtoll, final boolean dens) throws NRException {
        super(yy, dydxx, xx, atoll, rtoll, dens);
        dfdy = doub_mat(n, n);
        dfdx = doub_arr(n);
        k1 = doub_arr(n);
        k2 = doub_arr(n);
        k3 = doub_arr(n);
        k4 = doub_arr(n);
        k5 = doub_arr(n);
        k6 = doub_arr(n);
        cont1 = doub_arr(n);
        cont2 = doub_arr(n);
        cont3 = doub_arr(n);
        cont4 = doub_arr(n);
        a = doub_mat(n, n);
        // Input to the finalructor are the dependent variable y[0..n-1] and
        // its derivative dydx[0..n-1] at the starting value of the independent
        // variable x. Also input are the absolute and relative tolerances,
        // atol and rtol, and the boolean dense, which is true if dense output
        // is required.
        EPS = EPS(); // numeric_limits<double>::epsilon();
    }

    public void step(final double htry, final Dtype derivs) throws NRException {
        // Attempts a step with stepsize htry. On output, y and x are replaced
        // by their new values, hdid is the stepsize that was actually
        // accomplished, and hnext is the estimated next stepsize.
        double[] dydxnew = doub_arr(n);
        $double h = $(htry); // Set stepsize to the initial trial value.
        derivs.jacobian(x.$(), y.$(), dfdx, dfdy); // Compute the Jacobian and
                                                   // @f=@x.
        for (;;) {
            dy(h.$(), derivs); // Take a step.
            double err = error(); // Evaluate accuracy.
            if (con.success(err, h))
                break; // Step rejected. Try again with reduced h set
            if (abs(h.$()) <= abs(x.$()) * EPS) // by controller.
                throw new NRException("stepsize underflow in StepperRoss");
        }
        derivs.eval(x.$() + h.$(), yout, dydxnew); // Step succeeded.
        if (dense) // Compute coefficients for dense output.
            prepare_dense(h.$(), dydxnew);
        $$(dydx, dydxnew); // Reuse last derivative evaluation for next step.
        y.$(yout);
        xold = x.$(); // Used for dense output.
        x.$(x.$() + (hdid = h.$()));
        hnext = con.hnext;
    }

    public void dy(final double h, final Dtype derivs) throws NRException {
        // Given values for n variables y[0..n-1] and their derivatives
        // dydx[0..n-1] known at x, use the fourth-order stiffly stable
        // Rosenbrock method to advance the solution over an interval h and
        // store the incremented variables in yout[0..n-1]. Also store an
        // estimate of the local truncation error in yerr using the embedded
        // third-order method.
        double[] ytemp = doub_arr(n), dydxnew = doub_arr(n);
        int i;
        for (i = 0; i < n; i++) { // Set up the matrix 1=h  f 0.
            for (int j = 0; j < n; j++)
                a[i][j] = -dfdy[i][j];
            a[i][i] += 1.0 / (gam * h);
        }
        LUdcmp alu = new LUdcmp(a); // LU decomposition of the matrix.
        for (i = 0; i < n; i++)
            // Set up right-hand side for g1.
            ytemp[i] = dydx.$()[i] + h * d1 * dfdx[i];
        alu.solve(ytemp, k1); // Solve for g1.
        for (i = 0; i < n; i++)
            // Compute intermediate values of y.
            ytemp[i] = y.$()[i] + a21 * k1[i];
        derivs.eval(x.$() + c2 * h, ytemp, dydxnew); // Compute dydx at the
                                                     // intermediate values.
        for (i = 0; i < n; i++)
            // Set up right-hand side for g2.
            ytemp[i] = dydxnew[i] + h * d2 * dfdx[i] + c21 * k1[i] / h;
        alu.solve(ytemp, k2); // Solve for g2.
        for (i = 0; i < n; i++)
            // Compute intermediate values of y.
            ytemp[i] = y.$()[i] + a31 * k1[i] + a32 * k2[i];
        derivs.eval(x.$() + c3 * h, ytemp, dydxnew); // Compute dydx at the
                                                     // intermediate values.
        for (i = 0; i < n; i++)
            // Set up right-hand side for g3.
            ytemp[i] = dydxnew[i] + h * d3 * dfdx[i] + (c31 * k1[i] + c32 * k2[i]) / h;
        alu.solve(ytemp, k3); // Solve for g3.
        for (i = 0; i < n; i++)
            // Compute intermediate values of y.
            ytemp[i] = y.$()[i] + a41 * k1[i] + a42 * k2[i] + a43 * k3[i];
        derivs.eval(x.$() + c4 * h, ytemp, dydxnew); // Compute dydx at the
                                                     // intermediate values.
        for (i = 0; i < n; i++)
            // Set up right-hand side for g4.
            ytemp[i] = dydxnew[i] + h * d4 * dfdx[i] + (c41 * k1[i] + c42 * k2[i] + c43 * k3[i]) / h;
        alu.solve(ytemp, k4); // Solve for g4.
        for (i = 0; i < n; i++)
            // Compute intermediate values of y.
            ytemp[i] = y.$()[i] + a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i];
        double xph = x.$() + h;
        derivs.eval(xph, ytemp, dydxnew); // Compute dydx at the intermediate
                                          // values.
        for (i = 0; i < n; i++)
            // Set up right-hand side for g5.
            k6[i] = dydxnew[i] + (c51 * k1[i] + c52 * k2[i] + c53 * k3[i] + c54 * k4[i]) / h;
        alu.solve(k6, k5); // Solve for g5.
        for (i = 0; i < n; i++)
            // Compute the embedded solution.
            ytemp[i] += k5[i];
        derivs.eval(xph, ytemp, dydxnew); // Last derivative evaluation.
        for (i = 0; i < n; i++)
            // Compute the solution and the error.
            k6[i] = dydxnew[i] + (c61 * k1[i] + c62 * k2[i] + c63 * k3[i] + c64 * k4[i] + c65 * k5[i]) / h;
        alu.solve(k6, yerr);
        for (i = 0; i < n; i++)
            yout[i] = ytemp[i] + yerr[i];
    }

    public void prepare_dense(final double h, final double[] dydxnew) {
        // Store coefficients of interpolating polynomial for dense output
        // in cont1...cont4.
        for (int i = 0; i < n; i++) {
            cont1[i] = y.$()[i];
            cont2[i] = yout[i];
            cont3[i] = d21 * k1[i] + d22 * k2[i] + d23 * k3[i] + d24 * k4[i] + d25 * k5[i];
            cont4[i] = d31 * k1[i] + d32 * k2[i] + d33 * k3[i] + d34 * k4[i] + d35 * k5[i];
        }
    }

    public double dense_out(final int i, final double x, final double h) {
        // Evaluate interpolating polynomial for y[i] at location x, where
        // xold  x  xold Ch.
        double s = (x - xold) / h;
        double s1 = 1.0 - s;
        return cont1[i] * s1 + s * (cont2[i] + s1 * (cont3[i] + s * cont4[i]));
    }

    public double error() {

        // // Use yerr to compute norm of scaled error estimate. A value less
        // than one means the step was successful.
        double err = 0.0, sk;
        for (int i = 0; i < n; i++) {
            sk = atol + rtol * MAX(abs(y.$()[i]), abs(yout[i]));
            err += SQR(yerr[i] / sk);
        }
        return sqrt(err / n);
    }

    // Stepsize control depends on the fact that
    // yexact D y C O.h5/
    // yexact D yy C O.h4/
    // (17.5.34)
    // Thus jy  yyj D O.h4/ (17.5.35)
    // Referring back to the steps leading from equation (17.2.4) to equation
    // (17.2.12), we see that the new stepsize should be chosen as in
    // equation (17.2.12) but with the exponent 1/5 replaced by 1/4. Also,
    // experience shows that it is wise to prevent too large a stepsize change
    // in one step, otherwise we will probably have to undo the large change
    // in the next step. We adopt 0.2 and 6 as the maximum allowed decrease and
    // increase of h in one step. Methods for integrating stiff equations do
    // not suffer from the stability limitations that led to the PI controller
    // of 17.2.1. However, stiff problems often need a rapid decrease in
    // stepsize even when the previous step is successful. Also, sometimes the
    // effective order of the method can be lower than the simple Taylor series
    // prediction. Gustafsson [7] has proposed a predictive controller that
    // does a good job of dealing with these problems. The resulting formula
    // is
    // hnC1 D Shn
    // 
    // 1
    // errn
    // 1=4 hn
    // hn1
    // 
    // errn1
    // errn
    // 1=4
    // (17.5.36)
    // It is used only when a step is accepted.

    // struct Ross_constants { stepperross.h
    // Constants for the fourth-order stiy stable Rosenbrock method.
    // static const Doub c2,c3,c4,bet2p,bet3p,bet4p,d1,d2,d3,d4,a21,a31,a32,
    // a41,a42,a43,a51,a52,a53,a54,c21,c31,c32,c41,c42,c43,c51,c52,
    // c53,c54,c61,c62,c63,c64,c65,gam,d21,d22,d23,d24,d25,d31,d32,
    // d33,d34,d35;
    // };

    static final double c2 = 0.386;
    static final double c3 = 0.21;
    static final double c4 = 0.63;
    static final double bet2p = 0.0317;
    static final double bet3p = 0.0635;
    static final double bet4p = 0.3438;
    static final double d1 = 0.2500000000000000e+00;
    static final double d2 = -0.1043000000000000e+00;
    static final double d3 = 0.1035000000000000e+00;
    static final double d4 = -0.3620000000000023e-01;
    static final double a21 = 0.1544000000000000e+01;
    static final double a31 = 0.9466785280815826e+00;
    static final double a32 = 0.2557011698983284e+00;
    static final double a41 = 0.3314825187068521e+01;
    static final double a42 = 0.2896124015972201e+01;
    static final double a43 = 0.9986419139977817e+00;
    static final double a51 = 0.1221224509226641e+01;
    static final double a52 = 0.6019134481288629e+01;
    static final double a53 = 0.1253708332932087e+02;
    static final double a54 = -0.6878860361058950e+00;
    static final double c21 = -0.5668800000000000e+01;
    static final double c31 = -0.2430093356833875e+01;
    static final double c32 = -0.2063599157091915e+00;
    static final double c41 = -0.1073529058151375e+00;
    static final double c42 = -0.9594562251023355e+01;
    static final double c43 = -0.2047028614809616e+02;
    static final double c51 = 0.7496443313967647e+01;
    static final double c52 = -0.1024680431464352e+02;
    static final double c53 = -0.3399990352819905e+02;
    static final double c54 = 0.1170890893206160e+02;
    static final double c61 = 0.8083246795921522e+01;
    static final double c62 = -0.7981132988064893e+01;
    static final double c63 = -0.3152159432874371e+02;
    static final double c64 = 0.1631930543123136e+02;
    static final double c65 = -0.6058818238834054e+01;
    static final double gam = 0.2500000000000000e+00;
    static final double d21 = 0.1012623508344586e+02;
    static final double d22 = -0.7487995877610167e+01;
    static final double d23 = -0.3480091861555747e+02;
    static final double d24 = -0.7992771707568823e+01;
    static final double d25 = 0.1025137723295662e+01;
    static final double d31 = -0.6762803392801253e+00;
    static final double d32 = 0.6087714651680015e+01;
    static final double d33 = 0.1643084320892478e+02;
    static final double d34 = 0.2476722511418386e+02;
    static final double d35 = -0.6594389125716872e+01;
}
