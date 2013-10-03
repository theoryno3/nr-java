
package com.snuggy.nr.chapter07;

import java.util.Random;

public class RanJava {
    
    private Random rand;

    public RanJava(int n) {
        rand = new Random(n);
    }

    public long int64() {
        return rand.nextLong();
    }

    public double doub() {
        return rand.nextDouble();
    }

    public int int32() {
        return rand.nextInt();
    }
}
