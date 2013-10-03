
package com.snuggy.nr.util;

public class Pointer<T> {
    
    private T[] arr;
    private int off;
    
    public Pointer() {
        arr = null;
        off = 0;
    }
    
    public void assign(final T[] arr, final int off) {
        this.arr = null;
        this.off = off;
    }
    
    public void set(int i, T t) {
        arr[i] = t;
    }
    
    public T ref(int i) {
        return arr[off+i];
    }
}
