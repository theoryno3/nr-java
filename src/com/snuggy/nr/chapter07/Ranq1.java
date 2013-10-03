
package com.snuggy.nr.chapter07;

public class Ranq1 implements Random {

    // Recommended generator for everyday use. The period is
    // 1:8  1019. Calling conventions same as Ran, above.

    private long v;

    public Ranq1(long j) {
        v = 4101842887655102017L;
        v ^= j;
        v = int64();
    }

    public long int64() {
        v ^= v >>> 21;
        v ^= v << 35;
        v ^= v >>> 4;
        return v * 2685821657736338717L;
    }

    public double doub() {
        //return 5.42101086242752217E-20 * int64();
        long arg = int64();
        // Make sure it's in 0 to 1 rather than -0.5 to 0.5.
        // Save the right bit and shift to the right.
        int bit = (int) (arg & 0x0000000000000001L);
        arg >>>= 1;
        double r =  5.42101086242752217E-20 * arg;
        r *= 2.0;
        r += 5.42101086242752217E-20 * bit;
        return r;
    }

    public int int32() {
        return (int) int64();
    }

}
