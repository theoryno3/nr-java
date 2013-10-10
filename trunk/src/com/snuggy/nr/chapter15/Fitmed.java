
package com.snuggy.nr.chapter15;

import static com.snuggy.nr.chapter08.Static.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

public class Fitmed {

    // Object for fitting a straight line y D aCbx to a set of points .xi;
    // yi /, by the criterion of least absolute deviations. Call the constructor
    // to calculate the fit. The answers are then available as the variables
    // a, b, and abdev (the mean absolute deviation of the points from the
    // line).
    private int ndata;
    private double a, b, abdev; // Answers.
    private final double[] x, y;

    public double a() {
        return a;
    }

    public double b() {
        return b;
    }

    public double abdev() {
        return abdev;
    }

    public Fitmed(final double[] xx, final double[] yy) {
        // Constructor. Given a set of data points xx[0..ndata-1],
        // yy[0..ndata-1], sets a, b, and abdev.
        ndata = (xx.length);
        x = (xx);
        y = (yy);
        int j;
        double b1, b2, del, f, f1, f2 = 0.0, sigb, temp = 0.0;
        double sx = 0.0, sy = 0.0, sxy = 0.0, sxx = 0.0, chisq = 0.0;
        for (j = 0; j < ndata; j++) { // As a first guess for a and b, we will
                                      // find the
            sx += x[j]; // least-squares fitting line.
            sy += y[j];
            sxy += x[j] * y[j];
            sxx += SQR(x[j]);
        }
        del = ndata * sxx - sx * sx;
        a = (sxx * sy - sx * sxy) / del; // Least-squares solutions.
        b = (ndata * sxy - sx * sy) / del;
        for (j = 0; j < ndata; j++) {
            // chisq += (temp=y[j]-(a+b*x[j]), temp*temp); 
            // omfg are you kidding me? -- SKM
            temp = y[j] - (a + b * x[j]);
            chisq += (temp * temp);
        }
        sigb = sqrt(chisq / del); // The standard deviation will give some idea
                                  // of
        b1 = b; // how big an iteration step to take.
        f1 = rofunc(b1);
        if (sigb > 0.0) {
            b2 = b + SIGN(3.0 * sigb, f1); // Guess bracket as 3- away, in the
                                           // downhill di
            f2 = rofunc(b2); // rection known from f1.
            if (b2 == b1) {
                abdev /= ndata;
                return;
            }
            while (f1 * f2 > 0.0) { // Bracketing.
                b = b2 + 1.6 * (b2 - b1);
                b1 = b2;
                f1 = f2;
                b2 = b;
                f2 = rofunc(b2);
            }
            sigb = 0.01 * sigb;
            while (abs(b2 - b1) > sigb) {
                b = b1 + 0.5 * (b2 - b1); // Bisection.
                if (b == b1 || b == b2)
                    break;
                f = rofunc(b);
                if (f * f1 >= 0.0) {
                    f1 = f;
                    b1 = b;
                } else {
                    f2 = f;
                    b2 = b;
                }
            }
        }
        abdev /= ndata;
    }

    public double rofunc(final double b) {
        // Evaluates the right-hand side of equation (15.7.16) for a given
        // value of b.
        final double EPS = EPS(); // numeric_limits<Doub>::epsilon();
        int j;
        double d, sum = 0.0;
        final double[] arr = doub_vec(ndata);
        for (j = 0; j < ndata; j++)
            arr[j] = y[j] - b * x[j];
        if ((ndata & 1) == 1) {
            a = select((ndata - 1) >> 1, arr);
        } else {
            j = ndata >> 1;
            a = 0.5 * (select(j - 1, arr) + select(j, arr));
        }
        abdev = 0.0;
        for (j = 0; j < ndata; j++) {
            d = y[j] - (b * x[j] + a);
            abdev += abs(d);
            if (y[j] != 0.0)
                d /= abs(y[j]);
            if (abs(d) > EPS)
                sum += (d >= 0.0 ? x[j] : -x[j]);
        }
        return sum;
    }

}
