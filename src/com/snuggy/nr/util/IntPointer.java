
package com.snuggy.nr.util;

public class IntPointer {

    private final int[] arr;
    private final int off;

    public IntPointer(final int[] arr, final int off) {
        this.arr = arr;
        this.off = off;
    }

    public void set(int i, int t) {
        arr[off+i] = t;
    }

    public int ref(int i) {
        return arr[off+i];
    }

}
