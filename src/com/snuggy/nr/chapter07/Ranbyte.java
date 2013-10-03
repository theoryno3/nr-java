package com.snuggy.nr.chapter07;

import static com.snuggy.nr.util.Static.*;

public class Ranbyte {

    // Generator for random bytes using the algorithm generally known as RC4.
    private long[] s = long_arr(256);
    private int i, j;
    private long ss;
    private long v;

    public Ranbyte(long u) {
        // Constructor. Call with any integer seed.
        v = 2244614371L ^ u;
        for (i = 0; i < 256; i++) {
            s[i] = i;
        }
        for (j = 0, i = 0; i < 256; i++) {
            ss = s[i];
            j = (int) ((j + ss + (v >>> 24)) & 0xff);
            s[i] = s[j];
            s[j] = ss;
            v = (v << 24) | (v >>> 8);
        }
        i = (int) (j = 0);
        for (int k = 0; k < 256; k++)
            int8();
    }

    public int int8() {
        // Returns next random byte in the sequence.
        i = (i + 1) & 0xff;
        ss = s[i];
        j = (int) ((j + ss) & 0xff);
        s[i] = s[j];
        s[j] = ss;
        int r = (int) (s[(int) ((s[i] + s[j]) & 0xff)]);
        return r;
    }

    public int int32() {
        // Returns a random 32-bit integer constructed from 4 random bytes.
        // Slow!
        v = 0;
        for (int k = 0; k < 4; k++) {
            i = (i + 1) & 0xff;
            ss = s[i];
            j = (int) ((j + ss) & 0xff);
            s[i] = s[j];
            s[j] = ss;
            v = (v << 8) | s[(int) ((s[i] + s[j]) & 0xff)];
        }
        return (int) v;
    }

    public double doub() {
        // Returns a random double-precision floating value between 0. and 1.
        // Slow!!
        // return 2.32830643653869629E-10 *
        //          (int32() + 
        //           2.32830643653869629E-10 * int32());
        double r = 2.32830643653869629E-10 * 
                    ((int32() & 0xFFFFFFFFL) + 
                     2.32830643653869629E-10 * (int32() & 0xFFFFFFFFL));
        return r;
    }

}
