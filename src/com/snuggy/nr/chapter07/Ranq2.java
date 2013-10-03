
package com.snuggy.nr.chapter07;

public class Ranq2 implements Random {

    // Backup generator if Ranq1 has too short a period and Ran is too slow.
    // The period is 8:5  1037. Calling conventions same as Ran, above.

    private long v, w;

    public Ranq2(long j) {
        v = 4101842887655102017L;
        w = 1;
        v ^= j;
        w = int64();
        v = int64();
    }

    public long int64() {
        v ^= v >>> 17;
        v ^= v << 31;
        v ^= v >>> 8;
        w = 4294957665L*(w & 0xffffffffL) + (w >>> 32);
        return v ^ w;
    }

    public double doub() {
        //return 5.42101086242752217E-20 * int64();
        long arg = int64();
        // Make sure it's in 0 to 1 rather than -0.5 to 0.5.
        // Save the right bit and shift to the right.
        int bit = (int) (arg & 0x0000000000000001L);
        arg >>>= 1;
        double r = 5.42101086242752217E-20 * arg;
        r *= 2.0;
        r += 5.42101086242752217E-20 * bit;
        return r;
    }

    public int int32() {
        return (int) int64();
    }
}
