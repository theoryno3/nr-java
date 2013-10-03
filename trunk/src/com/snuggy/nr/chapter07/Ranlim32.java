
package com.snuggy.nr.chapter07;

public class Ranlim32 {

    // High-quality random generator using only 32-bit arithmetic. Same
    // conventions as Ran. Period 3:11  1037. Recommended only when
    // 64-bit arithmetic is not available.
    private int u;
    private int v = (int) (2244614371L);
    private int w1 = (521288629);
    private int w2 = (362436069);

    public Ranlim32(int j) {
        u = j ^ v;
        int32();
        v = u;
        int32();
    }

    public int int32() {
        u = (int) (u * 2891336453L + 1640531513);
        v ^= v >>> 13;
        v ^= v << 17;
        v ^= v >>> 5;
        w1 = 33378 * (w1 & 0xffff) + (w1 >>> 16);
        w2 = 57225 * (w2 & 0xffff) + (w2 >>> 16);
        int x = u ^ (u << 9);
        x ^= x >>> 17;
        x ^= x << 6;
        int y = w1 ^ (w1 << 17);
        y ^= y >>> 15;
        y ^= y << 5;
        return (int) ((x + v) ^ (y + w2));
    }

    public double doub() {
        //return 2.32830643653869629E-10 * int32();
        int arg = int32();
        double r = scale(arg);
        return r;
    }

    public double truedoub() {
        double t = scale(int32());
        return scale(int32()) +
               2.32830643653869629E-10 * t;
    }
    
    private static double scale(int arg) {
        // Make sure it's in 0 to 1 rather than -0.5 to 0.5.
        // Save the right bit and shift to the right.
        int bit = (arg & 0x00000001);
        arg >>>= 1;
        double r = 2.32830643653869629E-10 * arg;
        r *= 2.0;
        r += 2.32830643653869629E-10 * bit;
        return r;
    }

}
