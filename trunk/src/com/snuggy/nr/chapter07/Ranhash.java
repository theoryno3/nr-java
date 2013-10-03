
package com.snuggy.nr.chapter07;

import java.nio.*;

public class Ranhash {

    // High-quality random hash of an integer into several numeric types.

    public long int64(final LongBuffer u_arr, final int u_off) {
        return int64(u_arr.get(u_off));
    }

    public long int64(final long[] u_arr, final int u_off) {
        return int64(u_arr[u_off]);
    }

    public long int64(final long u) {
        // Returns hash of u as a 64-bit integer.
        long v = u * 3935559000370003845L + 2691343689449507681L;
        v ^= v >>> 21;
        v ^= v << 37;
        v ^= v >>> 4;
        v *= 4768777513237032717L;
        v ^= v << 20;
        v ^= v >>> 41;
        v ^= v << 5;
        return v;
    }

    public int int32(final LongBuffer u_arr, final int u_off) {
        return int32(u_arr.get(u_off));
    }

    public int int32(final long[] u_arr, final int u_off) {
        return int32(u_arr[u_off]);
    }

    public int int32(final long u) {
        // Returns hash of u as a 32-bit integer.
        return (int) (int64(u) & 0xffffffffL);
    }

    public double doub(final LongBuffer u_arr, final int u_off) {
        return doub(u_arr.get(u_off));
    }

    public double doub(final long[] u_arr, final int u_off) {
        return doub(u_arr[u_off]);
    }

    public double doub(final long u) {
        // Returns hash of u as a double-precision floating value 
        // between 0. and 1.
        //return 5.42101086242752217E-20 * int64(u);
        // Make sure it's in 0 to 1 rather than -0.5 to 0.5.
        long arg = int64(u);
        int bit = (int) (arg & 0x0000000000000001L);
        arg >>>= 1;
        double r =  5.42101086242752217E-20 * arg;
        r *= 2.0;
        r += 5.42101086242752217E-20 * bit;
        return r;
    }

}
