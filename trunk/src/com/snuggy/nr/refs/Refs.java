
package com.snuggy.nr.refs;

import com.snuggy.nr.util.*;

public class Refs {
    
    // references to ints
    
    public static $int $(final int t) {
        $int r = new IntRef(t);
        return r;
    }
    
    public static int $(final $int x, int y) {
        x.$(y);
        return y;
    }
    
    public static int $(final $int x, final $int y) {
        x.$(y.$());
        return y.$();
    }
    
    // references to doubles
    
    public static $double $(final double t) {
        $double r = new DoubleRef(t);
        return r;
    }
    
    public static $double $(final double[] arr, final int off) {
        $double r = new DoubleRefFromArrayElement(arr, off);
        return r;
    }
    
    public static double $(final $double x, double y) {
        x.$(y);
        return y;
    }
    
    public static double $(final $double x, final $double y) {
        x.$(y.$());
        return y.$();
    }
    
    // references to booleans
    
    public static $boolean $(final boolean t) {
        $boolean r = new BooleanRef(t);
        return r;
    }
    
    public static void $(final $boolean x, boolean y) {
        x.$(y);
    }
    
    public static void $(final $boolean x, final $boolean y) {
        x.$(y.$());
    }
    
    // references to parameterized Objects passed by reference
    
	public static <T> $<T> $(final T t) {
	    $<T> r = new ObjectRef<T>(t);
	    return r;
	}
	
	public static <T> void $(final $<T> x, final T y) {
	    x.$(y);
	}
	
	public static <T> void $(final $<T> x, final $<T> y) {
	    x.$(y.$());
	}
	
    // references to parameterized Objects passed by reference or value
    
	public static <T extends ByValue<T>> $$<T> $$(final T t) {
	    $$<T> r = new ObjectRefByValue<T>(t);
	    return r;
	}
	
	public static <T extends ByValue<T>> void $$(final $$<T> x, final T y) {
	    x.$$(y);
	}
	
	public static <T extends ByValue<T>> void $$(final $$<T> x, final $<T> y) {
	    x.$$(y.$());
	}
	
	public static <T extends ByValue<T>> void $$(final $<T> x, final $$<T> y) {
	    x.$(y.$$());
	}
	
    // references to double[] 
    
	public static $double1d $(final double[] t) {
	    $double1d r = new Double1DWrapper(t);
	    return r;
	}
	
	public static $double1d $$(final double[] t) throws NRException {
	    $double1d r = $(new Double1DWrapper(t).$$());
	    return r;
	}
	
	public static void $$(final $double1d x, final double[] y) throws NRException {
	    x.$$(y);
	}
	
	public static void $$(final $double1d x, final $double1d y) throws NRException {
	    x.$$(y.$());
	}
	
    // references to double[][] 
    
	public static $double2d $(final double[][] t) {
	    $double2d r = new Double2DWrapper(t);
	    return r;
	}
	
	public static $double2d $$(final double[][] t) {
	    $double2d r = $(new Double2DWrapper(t).$$());
	    return r;
	}
	
	public static void $$(final $double2d x, final $double2d y) throws NRException {
	    x.$$(y.$());
	}
	
	public static void $$(final $double2d x, final double[][] y) throws NRException {
	    x.$$(y);
	}
	
    // references to int[] 
    
	public static $int1d $(final int[] t) {
	    $int1d r = new Int1DWrapper(t);
	    return r;
	}
	
	public static $int1d $$(final int[] t) {
	    $int1d r = $(new Int1DWrapper(t).$$());
	    return r;
	}
	
	public static void $$(final $int1d x, final $int1d y) throws NRException {
	    x.$$(y.$());
	}
	
	public static void $$(final $int1d x, final int[] y) throws NRException {
	    x.$$(y);
	}
	
    // references to int[][] 
    
	public static $int2d $(final int[][] t) {
	    $int2d r = new Int2DWrapper(t);
	    return r;
	}
	
	public static $int2d $$(final int[][] t) {
	    $int2d r = $(new Int2DWrapper(t).$$());
	    return r;
	}
	
	public static void $$(final $int2d x, final int[][] y) throws NRException {
	    x.$$(y);
	}
	
	public static void $$(final $int2d x, final $int2d y) throws NRException {
	    x.$$(y.$());
	}
	
	// classes
	
    static class Double1DWrapper implements $double1d {
        private double[] t;
        public Double1DWrapper(final double[] t) {
            this.t = t;
        }
        @Override
        public void $(final double[] t) {
            this.t = t;
        }
        @Override
        public double[] $() {
            return t;
        }
        @Override
        public double[] $$() {
            double[] r = new double[t.length];
            System.arraycopy(t, 0, r, 0, t.length);
            return r;
        }
        @Override
        public void $$(double[] t) throws NRException {
            if (this.t.length != t.length)
                throw new NRException("this.t.length != t.length");
            System.arraycopy(t, 0, this.t, 0, t.length);
        }
        @Override
        public String toString() {
            return t.toString();
        }
    }
	
    static class Double2DWrapper implements $double2d {
        private double[][] t;
        public Double2DWrapper(final double[][] t) {
            this.t = t;
        }
        @Override
        public void $(final double[][] t) {
            this.t = t;
        }
        @Override
        public double[][] $() {
            return t;
        }
        @Override
        public double[][] $$() {
            double[][] r = new double[t.length][];
            for (int i = 0; i < t.length; i++) {
                r[i] = new double[t[i].length];
	            System.arraycopy(t[i], 0, r[i], 0, t[i].length);
            }
            return r;
        }
        @Override
        public void $$(double[][] t) throws NRException {
            if (this.t.length != t.length)
                throw new NRException("this.t.length != t.length");
            for (int i = 0; i < t.length; i++) {
                if (this.t[i].length != t[i].length)
                    throw new NRException("this.t[i].length !+ t[i].length");
	            System.arraycopy(t[i], 0, this.t[i], 0, t[i].length);
            }
        }
        @Override
        public String toString() {
            return t.toString();
        }
    }
	
    static class Int1DWrapper implements $int1d {
        private int[] t;
        public Int1DWrapper(final int[] t) {
            this.t = t;
        }
        @Override
        public void $(final int[] t) {
            this.t = t;
        }
        @Override
        public int[] $() {
            return t;
        }
        @Override
        public int[] $$() {
            int[] r = new int[t.length];
            System.arraycopy(t, 0, r, 0, t.length);
            return r;
        }
        @Override
        public void $$(int[] t) throws NRException {
            if (this.t.length != t.length)
                throw new NRException("this.t.length != t.length");
            System.arraycopy(t, 0, this.t, 0, t.length);
        }
        @Override
        public String toString() {
            return t.toString();
        }
    }
	
    static class Int2DWrapper implements $int2d {
        private int[][] t;
        public Int2DWrapper(final int[][] t) {
            this.t = t;
        }
        @Override
        public void $(final int[][] t) {
            this.t = t;
        }
        @Override
        public int[][] $() {
            return t;
        }
        @Override
        public int[][] $$() {
            int[][] r = new int[t.length][];
            for (int i = 0; i < t.length; i++) {
                r[i] = new int[t[i].length];
	            System.arraycopy(t[i], 0, r[i], 0, t[i].length);
            }
            return r;
        }
        @Override
        public void $$(int[][] t) throws NRException {
            if (this.t.length != t.length)
                throw new NRException("this.t.length != t.length");
            for (int i = 0; i < t.length; i++) {
                if (this.t[i].length != t[i].length)
                    throw new NRException("this.t[i].length !+ t[i].length");
	            System.arraycopy(t[i], 0, this.t[i], 0, t[i].length);
            }
        }
        @Override
        public String toString() {
            return t.toString();
        }
    }
	
	static class DoubleRef implements $double {
        private double t;
        public DoubleRef(double t) {
            this.t = t;
        }
        @Override
        public void $(double t) {
            this.t = t;
        }
        @Override
        public double $() {
            return t;
        }
        @Override
        public String toString() {
            return String.valueOf(t);
        }
	}
	
	static class DoubleRefFromArrayElement implements $double {
        private double[] arr;
        private int off;
        public DoubleRefFromArrayElement(final double[] arr, final int off) {
            this.arr = arr;
            this.off = off;
        }
        @Override
        public void $(double t) {
            arr[off] = t;
        }
        @Override
        public double $() {
            return arr[off];
        }
        @Override
        public String toString() {
            return String.valueOf(arr[off]);
        }
	}
	
	static class IntRef implements $int {
        private int t;
        public IntRef(int t) {
            this.t = t;
        }
        @Override
        public void $(int t) {
            this.t = t;
        }
        @Override
        public int $() {
            return t;
        }
        @Override
        public String toString() {
            return String.valueOf(t);
        }
	}
	
	static class BooleanRef implements $boolean {
        private boolean t;
        public BooleanRef(boolean t) {
            this.t = t;
        }
        @Override
        public void $(boolean t) {
            this.t = t;
        }
        @Override
        public boolean $() {
            return t;
        }
        @Override
        public String toString() {
            return String.valueOf(t);
        }
	}
	
    static class ObjectRef<T> implements $<T> {
        private T t;
        public ObjectRef(T t) {
            this.t = t;
        }
        @Override
        public void $(T t) {
            this.t = t;
        }
        @Override
        public T $() {
            return t;
        }
        @Override
        public String toString() {
            return t.toString();
        }
    }
	
    static class ObjectRefByValue<T extends ByValue<T>> implements $$<T> {
        private T t;
        public ObjectRefByValue(T t) {
            this.t = t;
        }
        @Override
        public T $() {
            return t;
        }
        @Override
        public void $(T t) {
            this.t = t;
        }
        @Override
        public T $$() {
            return t.copyOut();
        }
        @Override
        public void $$(T t) {
            this.t.copyIn(t);
        }
        @Override
        public String toString() {
            return t.toString();
        }
    }
}
