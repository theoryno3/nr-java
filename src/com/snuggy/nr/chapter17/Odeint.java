
package com.snuggy.nr.chapter17;

import static com.snuggy.nr.refs.Refs.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import java.lang.reflect.*;

import com.snuggy.nr.refs.*;
import com.snuggy.nr.util.*;

@Deprecated @Broken
public class Odeint {

    // Driver for ODE solvers with adaptive stepsize control. The template
    // parameter should be one of the derived classes of StepperBase defining
    // a particular integration algorithm.
    private static final int MAXSTP = 50000; // Take at most MAXSTP steps.
    private double EPS;
    @SuppressWarnings("unused")
    private int nok;
    @SuppressWarnings("unused")
    private int nbad;
    private int nvar;
    private double x1, x2, hmin;
    private boolean dense; // true if dense output requested by
    private final double[] y, dydx; // out.
    private final double[] ystart;
    private final Output out;

    private Dtype derivs; // Get the type of derivs from the
    private IStepperBS s; // stepper.
    private int nstp;
    private $double x;
    private double h;

    // public Odeint(final double[] ystartt, final double xx1, final double xx2,
    // final double atol, final double rtol, final double h1,
    // final double hminn, final Output outt, final HasDenseOut derivss);
    // Constructor sets everything up. The routine integrates starting values
    // ystart[0..nvar-1] from xx1 to xx2 with absolute tolerance atol and
    // relative tolerance rtol. The quantity h1 should be set as a guessed
    // first stepsize, hmin as the minimum allowed stepsize (can be zero).
    // An Output object out should be input to control the saving of
    // intermediate values.
    // On output, nok and nbad are the number of good and bad (but retried and
    // fixed) steps taken, and ystart is replaced by values at the end of the
    // integration interval. derivs is the user-supplied routine (function or
    // functor) for calculating the right-hand side derivative.
    // void integrate(); Does the actual integration.

    public Odeint(final Class<? extends StepperBS> className, final double[] ystartt, final double xx1,
            final double xx2, final double atol, final double rtol, final double h1, final double hminn,
            final Output outt, final Dtype derivss) throws InstantiationException, IllegalAccessException, NRException,
            NoSuchMethodException, SecurityException, IllegalArgumentException, InvocationTargetException {
        nvar = (ystartt.length);
        y = doub_arr(nvar);
        dydx = doub_arr(nvar);
        ystart = (ystartt);
        x = $(xx1);
        nok = (0);
        nbad = (0);
        x1 = (xx1);
        x2 = (xx2);
        hmin = (hminn);
        dense = (outt.dense());
        out = (outt);
        // derivs = (derivss);
        // s = (y,dydx,x,atol,rtol,dense);
        Constructor<? extends IStepperBS> constructor = className.getConstructor(double[].class, double[].class,
                $double.class, double.class, double.class, boolean.class);
        s = constructor.newInstance(y, dydx, x, atol, rtol, dense);
        derivs = derivss;
        EPS = EPS(); // numeric_limits<Doub>::epsilon();
        h = SIGN(h1, x2 - x1);
        for (int i = 0; i < nvar; i++)
            y[i] = ystart[i];
        out.init(s.neqn(), x1, x2);
    }

    public void integrate() throws NRException {
        derivs.eval(x.$(), y, dydx);
        if (dense) // Store initial values.
            out.out(-1, x.$(), y, s, h);
        else
            out.save(x.$(), y);
        for (nstp = 0; nstp < MAXSTP; nstp++) {
            if ((x.$() + h * 1.0001 - x2) * (x2 - x1) > 0.0)
                h = x2 - x.$(); // If stepsize can overshoot, decrease.
            s.step(h, derivs); // Take a step.
            if (s.hdid() == h)
                ++nok;
            else
                ++nbad;
            if (dense)
                out.out(nstp, x.$(), y, s, s.hdid());
            else
                out.save(x.$(), y);
            if ((x.$() - x2) * (x2 - x1) >= 0.0) { // Are we done?
                for (int i = 0; i < nvar; i++)
                    ystart[i] = y[i]; // Update ystart.
                if (out.kmax() > 0 && abs(out.xsave()[out.count() - 1] - x2) > 100.0 * abs(x2) * EPS)
                    out.save(x.$(), y); // Make sure last step gets saved.
                return; // Normal exit.
            }
            if (abs(s.hnext()) <= hmin)
                throw new NRException("Step size too small in Odeint");
            h = s.hnext();
        }
        throw new NRException("Too many steps in routine Odeint");
    }

}
