
package com.snuggy.nr.chapter05;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.util.*;

public class Eulsum {

    // Convergence acceleration of an alternating series by the Euler
    // transformation. Initialize by calling the constructor with arguments
    // nmax, an upper bound on the number of terms to be summed, and epss,
    // the desired accuracy. Then make successive calls to the function next
    // (see below).
    private final double[] wksp;
    private int n, ncv;
    private boolean cnvgd;
    private double sum, eps, lastval, lasteps;

    public Eulsum(final int nmax, final double epss) {
        wksp = doub_vec(nmax);
        n = (0);
        ncv = (0);
        cnvgd = (false);
        sum = (0.);
        eps = (epss);
        lastval = (0.);
    }
    
    public boolean cnvgd() {
        return cnvgd;
    }

    public double next(final double term) throws NRException {
        // Incorporates into sum the next term, with value term, of an
        // alternating series. On each call term should have a sign opposite
        // to that of the previous call. The ag cnvgd is set when convergence
        // is detected.
        int j;
        double tmp, dum;
        if (n + 1 > wksp.length)
            throw new NRException("wksp too small in eulsum");
        if (n == 0) { // Initialize:
            sum = 0.5 * (wksp[n++] = term); // Return rst estimate.
        } else {
            tmp = wksp[0];
            wksp[0] = term;
            for (j = 1; j < n; j++) { // Update saved quantities by van Wijn
                dum = wksp[j]; // gaarden's algorithm.
                wksp[j] = 0.5 * (wksp[j - 1] + tmp);
                tmp = dum;
            }
            wksp[n] = 0.5 * (wksp[n - 1] + tmp);
            if (abs(wksp[n]) <= abs(wksp[n - 1])) // Favorable to increase p,
                sum += (0.5 * wksp[n++]); // and the table becomes longer.
            // else Favorable to increase n,
            sum += wksp[n]; // the table doesn't become longer.
        }
        lasteps = abs(sum - lastval);
        if (lasteps <= eps)
            ncv++;
        if (ncv >= 2)
            cnvgd = true;
        return (lastval = sum);
    }

}
