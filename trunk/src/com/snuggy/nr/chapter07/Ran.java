
package com.snuggy.nr.chapter07;

public class Ran implements Random {

    // Implementation of the highest quality recommended generator. The
    // constructor is called with an integer seed and creates an instance of
    // the generator. The member functions int64, doub, and int32 return the
    // next values in the random sequence, as a variable type indicated by
    // their names. The period of the generator is 3:138  1057.

    private long u;
    private long v = 4101842887655102017L;
    private long w = 1;

    public Ran(final long j) {
        // Constructor. Call with any integer seed (except value of v above).
        u = j ^ v;
        int64();
        v = u;
        int64();
        w = v;
        int64();
    }

    public long int64() {
        u = u * 2862933555777941757L + 7046029254386353087L;
        v ^= v >>> 17;
        v ^= v << 31;
        v ^= v >>> 8;
        w = 4294957665L * (w & 0xffffffffL) + (w >>> 32);
        long x = u ^ (u << 21);
        x ^= x >>> 35;
        x ^= x << 4;
        return (x + v) ^ w;
    }

    public double doub() {
        // return 5.42101086242752217E-20 * int64();
        /*
        double r = 5.42101086242752217E-20 * int64();
        if (r < 0.0)
            r += 1.0;
        return r;
        */
        long arg = int64();
        // Make sure it's in 0 to 1 rather than -0.5 to 0.5.
        // Save the right bit and shift to the right.
        int bit = (int) (arg & 0x0000000000000001L);
        arg >>>= 1;
        double r =  5.42101086242752217E-20 * arg;
        r *= 2.0;
        r += 5.42101086242752217E-20 * bit;
        return r;
    } // Return random double-precision floating value in the range 0. to 1.

    public int int32() {
        int r = (int) int64();
        return r;
    }
    // Return 32-bit random integer.

}
