
package com.snuggy.nr.chapter07;

import static com.snuggy.nr.util.Static.*;

public class Ranfib {

    // Implements Knuth’s subtractive generator using only floating operations.
    // See text for cautions.
    private double[] dtab = doub_arr(55);
    private double dd;
    private int inext, inextp;

    public Ranfib(long j) {
        // Constructor. Call with any integer seed. Uses Ranq1 to initialize.
        inext = (0);
        inextp = (31);
        Ranq1 init = new Ranq1(j);
        for (int k = 0; k < 55; k++)
            dtab[k] = init.doub();
    }

    public double doub() {
        // Returns random double-precision floating value between 0. and 1.
        if (++inext == 55)
            inext = 0;
        if (++inextp == 55)
            inextp = 0;
        dd = dtab[inext] - dtab[inextp];
        if (dd < 0)
            dd += 1.0;
        return (dtab[inext] = dd);
    }

    public long int32() {
        // Returns random 32-bit integer. Recommended only for testing purposes.
        return (long) (doub() * 4294967295.0);
    }

}
