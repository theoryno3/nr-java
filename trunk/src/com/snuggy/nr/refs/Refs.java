
package com.snuggy.nr.refs;

import com.snuggy.nr.util.*;

public class Refs {
    
    // reference to int
    
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
    
    // reference to int[]
    
    public static $int1d $(final int[] t) {
        $int1d r = new Int1DRef(t);
        return r;
    }
    
    public static void $(final $int1d x, int[] y) {
        x.$(y);
    }
    
    public static void $(final $int1d x, final $int1d y) {
        x.$(y.$());
    }
    
    // reference to int[][]
    
    public static $int2d $(final int[][] t) {
        $int2d r = new Int2DRef(t);
        return r;
    }
    
    public static void $(final $int2d x, int[][] y) {
        x.$(y);
    }
    
    public static void $(final $int2d x, final $int2d y) {
        x.$(y.$());
    }
    
    // reference to double
    
    public static $double $(final double t) {
        $double r = new DoubleRef(t);
        return r;
    }
    
    public static $double $(final double[] arr, final int off) {
        $double r = new DoubleRefFromArrayElement(arr, off);
        return r;
    }
    
    public static $double $(final $double x, final int off) throws NRException {
        $double r = new DoubleRefFromArrayElement(x, off);
        return r;
    }
    
    public static $double $(final double[][] arr, final int i, final int j) {
        $double r = new DoubleRefFrom2DArrayElement(arr, i, j);
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
    
    // reference to boolean
    
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
    
    // reference to double[]
    
    public static $double1d $(final double[] t) {
        $double1d r = new Double1DRef(t);
        return r;
    }
    
    public static void $(final $double1d x, double[] y) {
        x.$(y);
    }
    
    public static void $(final $double1d x, final $double1d y) {
        x.$(y.$());
    }
    
    // reference to double[][]
    
    public static $double2d $(final double[][] t) {
        $double2d r = new Double2DRef(t);
        return r;
    }
    
    public static void $(final $double2d x, double[][] y) {
        x.$(y);
    }
    
    public static void $(final $double2d x, final $double2d y) {
        x.$(y.$());
    }
    
    // reference to parameterized Object passed by reference
    
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
	
	public static <T extends ByValue<T>> void $$(final $$<T> x, final T y) throws NRException {
	    x.$$(y);
	}
	
	public static <T extends ByValue<T>> void $$(final $$<T> x, final $<T> y) throws NRException {
	    x.$$(y.$());
	}
	
	public static <T extends ByValue<T>> void $$(final $<T> x, final $$<T> y) throws NRException {
	    x.$(y.$$());
	}
	
    // double[] by value
    
	public static $$double1d $$(final double[] t) throws NRException {
	    $$double1d r = new Double1DByVal(new double[t.length]);
	    $$(r, t);
	    return r;
	}
	
	/*
	public static $$double1d $$(final double[] arr[], final int off) {
	    $$double1d r = new Double1DRefFromArrayElement(arr, off);
	    return r;
	}
	*/
	
	public static void $$(final $$double1d x, final double[] y) throws NRException {
	    x.$$(y);
	}
	
	public static void $$(final $$double1d x, final $$double1d y) throws NRException {
	    x.$$(y.$());
	}
	
	public static void $$(final double[] x, final $$double1d y) throws NRException {
	    y.copyOut(x);
	}
	
	public static void $$(final double[] x, final double[] y) throws NRException {
	    if (x.length != y.length)
	        throw new NRException("x.length != y.length");
	    System.arraycopy(y, 0, x, 0, x.length);
	}
	
	public static void $(final $$double1d x, final int i, final double v) throws NRException {
	    x.$(i, v);
	}
	
    // double[][] by value
    
	public static $$double2d $$(final double[][] t) throws NRException {
	    $$double2d r = new Double2DByVal(t.length, t[0].length);
	    $$(r, t);
	    return r;
	}
	
	public static void $$(final $$double2d x, final $$double2d y) throws NRException {
	    x.$$(y.$());
	}
	
	public static void $$(final $$double2d x, final double[][] y) throws NRException {
	    x.$$(y);
	}
	
	public static double $(final $$double2d x, final int i, final int j, double v) {
	    x.$(i, j, v);
	    return v;
	}
	
    // int[] by value
    
	public static $$int1d $$(final int[] t) {
	    $$int1d r = new Int1DByVal(t);
	    return r;
	}
	
	public static void $$(final $$int1d x, final $$int1d y) throws NRException {
	    x.$$(y.$());
	}
	
	public static void $$(final $$int1d x, final int[] y) throws NRException {
	    x.$$(y);
	}
	
    // int[][] by value
    
	public static $$int2d $$(final int[][] t) {
	    $$int2d r = new Int2DByVal(t);
	    return r;
	}
	
	public static void $$(final $$int2d x, final int[][] y) throws NRException {
	    x.$$(y);
	}
	
	public static void $$(final $$int2d x, final $$int2d y) throws NRException {
	    x.$$(y.$());
	}
	
    // boolean[] by value
    
	public static $$boolean1d $$(final boolean[] t) {
	    $$boolean1d r = new Boolean1DByVal(t);
	    return r;
	}
	
	public static void $$(final $$boolean1d x, final $$boolean1d y) throws NRException {
	    x.$$(y.$());
	}
	
	public static void $$(final $$boolean1d x, final boolean[] y) throws NRException {
	    x.$$(y);
	}
	
    // boolean[][] by value
    
	public static $$boolean2d $$(final boolean[][] t) {
	    $$boolean2d r = new Boolean2DByVal(t);
	    return r;
	}
	
	public static void $$(final $$boolean2d x, final boolean[][] y) throws NRException {
	    x.$$(y);
	}
	
	public static void $$(final $$boolean2d x, final $$boolean2d y) throws NRException {
	    x.$$(y.$());
	}
	
	// classes
	
    static class Double1DByVal implements $$double1d {
        private double[] t;
        public Double1DByVal(final double[] t) {
            this.t = t;
        }
        @Override
        public void $(final double[] t) {
            this.t = t;
        }
        @Override
        public final double[] $() {
            return t;
        }
        @Override
        public final double[] $$() {
            final double[] r = new double[t.length];
            System.arraycopy(t, 0, r, 0, t.length);
            return r;
        }
        @Override
        public void $$(final double[] t) throws NRException {
            if (this.t.length != t.length)
                throw new NRException("this.t.length != t.length");
            System.arraycopy(t, 0, this.t, 0, t.length);
        }
        @Override
        public String toString() {
            return t.toString();
        }
        @Override
        public void $(int i, double v) {
            t[i] = v;
        }
        @Override
        public void copyOut(double[] t) throws NRException {
            if (this.t.length != t.length)
                throw new NRException("this.t.length != t.length");
            System.arraycopy(this.t, 0, t, 0, t.length);
        }
    }
	
    static class Double1DByValueFromArrayElement implements $$double1d {
        private final double[] arr[];
        private final int off;
        public Double1DByValueFromArrayElement(final double[] arr[], final int off) {
            this.arr = arr;
            this.off = off;
        }
        @Override
        public void $(final double[] t) {
            arr[off] = t;
        }
        @Override
        public final double[] $() {
            return arr[off];
        }
        @Override
        public final double[] $$() {
            final double[] r = new double[arr[off].length];
            System.arraycopy(arr[off], 0, r, 0, arr[off].length);
            return r;
        }
        @Override
        public void $$(final double[] t) throws NRException {
            if (this.arr[off].length != arr[off].length)
                throw new NRException("this.t.length != arr[off].length");
            System.arraycopy(t, 0, this.arr[off], 0, arr[off].length);
        }
        @Override
        public String toString() {
            return arr[off].toString();
        }
        @Override
        public void $(int i, double v) throws NRException {
            throw new NRException();
        }
        @Override
        public void copyOut(double[] t) throws NRException {
            throw new NRException();
        }
    }
	
    static class Double2DByVal implements $$double2d {
        private double[][] t;
        public Double2DByVal(final int m, final int n) {
            this.t = new double[m][n];
            for (int i = 0; i < t.length; i++) {
	            System.arraycopy(t[i], 0, this.t[i], 0, t[i].length);
            }
        }
        @Override
        public void $(final double[][] t) {
            this.t = t;
        }
        @Override
        public final double[][] $() {
            return t;
        }
        @Override
        public final double[][] $$() {
            final double[][] r = new double[t.length][];
            for (int i = 0; i < t.length; i++) {
                r[i] = new double[t[i].length];
	            System.arraycopy(t[i], 0, r[i], 0, t[i].length);
            }
            return r;
        }
        @Override
        public void $$(final double[][] t) throws NRException {
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
        @Override
        public void $(int i, int j, double v) {
            t[i][j] = v;
        }
    }
	
    static class Int1DByVal implements $$int1d {
        private int[] t;
        public Int1DByVal(final int[] t) {
            this.t = t;
        }
        @Override
        public void $(final int[] t) {
            this.t = t;
        }
        @Override
        public final int[] $() {
            return t;
        }
        @Override
        public final int[] $$() {
            final int[] r = new int[t.length];
            System.arraycopy(t, 0, r, 0, t.length);
            return r;
        }
        @Override
        public void $$(final int[] t) throws NRException {
            if (this.t.length != t.length)
                throw new NRException("this.t.length != t.length");
            System.arraycopy(t, 0, this.t, 0, t.length);
        }
        @Override
        public String toString() {
            return t.toString();
        }
    }
	
    static class Int2DByVal implements $$int2d {
        private int[][] t;
        public Int2DByVal(final int[][] t) {
            this.t = t;
        }
        @Override
        public void $(final int[][] t) {
            this.t = t;
        }
        @Override
        public final int[][] $() {
            return t;
        }
        @Override
        public final int[][] $$() {
            final int[][] r = new int[t.length][];
            for (int i = 0; i < t.length; i++) {
                r[i] = new int[t[i].length];
	            System.arraycopy(t[i], 0, r[i], 0, t[i].length);
            }
            return r;
        }
        @Override
        public void $$(final int[][] t) throws NRException {
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
	
    static class Boolean1DByVal implements $$boolean1d {
        private boolean[] t;
        public Boolean1DByVal(final boolean[] t) {
            this.t = t;
        }
        @Override
        public void $(final boolean[] t) {
            this.t = t;
        }
        @Override
        public final boolean[] $() {
            return t;
        }
        @Override
        public final boolean[] $$() {
            final boolean[] r = new boolean[t.length];
            System.arraycopy(t, 0, r, 0, t.length);
            return r;
        }
        @Override
        public void $$(final boolean[] t) throws NRException {
            if (this.t.length != t.length)
                throw new NRException("this.t.length != t.length");
            System.arraycopy(t, 0, this.t, 0, t.length);
        }
        @Override
        public String toString() {
            return t.toString();
        }
    }
	
    static class Boolean2DByVal implements $$boolean2d {
        private boolean[][] t;
        public Boolean2DByVal(final boolean[][] t) {
            this.t = t;
        }
        @Override
        public void $(final boolean[][] t) {
            this.t = t;
        }
        @Override
        public final boolean[][] $() {
            return t;
        }
        @Override
        public final boolean[][] $$() {
            final boolean[][] r = new boolean[t.length][];
            for (int i = 0; i < t.length; i++) {
                r[i] = new boolean[t[i].length];
	            System.arraycopy(t[i], 0, r[i], 0, t[i].length);
            }
            return r;
        }
        @Override
        public void $$(final boolean[][] t) throws NRException {
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
        @Override
        public double $_(int n) throws NRException {
            throw new NRException("offset used for non-offset $double");
        }
        @Override
        public void $_(int n, double v) throws NRException {
            throw new NRException("offset used for non-offset $double");
        }
	}
	
	static class Double1DRef implements $double1d {
        private double[] t;
        public Double1DRef(double[] t) {
            this.t = t;
        }
        @Override
        public void $(double[] t) {
            this.t = t;
        }
        @Override
        public double[] $() {
            return t;
        }
        @Override
        public String toString() {
            return String.valueOf(t);
        }
	}
	
	static class Int1DRef implements $int1d {
        private int[] t;
        public Int1DRef(int[] t) {
            this.t = t;
        }
        @Override
        public void $(int[] t) {
            this.t = t;
        }
        @Override
        public int[] $() {
            return t;
        }
        @Override
        public String toString() {
            return String.valueOf(t);
        }
	}
	
	static class Int2DRef implements $int2d {
        private int[][] t;
        public Int2DRef(int[][] t) {
            this.t = t;
        }
        @Override
        public void $(int[][] t) {
            this.t = t;
        }
        @Override
        public int[][] $() {
            return t;
        }
        @Override
        public String toString() {
            return String.valueOf(t);
        }
	}
	
	static class Double2DRef implements $double2d {
        private double[][] t;
        public Double2DRef(double[][] t) {
            this.t = t;
        }
        @Override
        public void $(double[][] t) {
            this.t = t;
        }
        @Override
        public double[][] $() {
            return t;
        }
        @Override
        public String toString() {
            return String.valueOf(t);
        }
	}
	
	static class DoubleRefFromArrayElement implements $double {
        private final double[] arr;
        private int off;
        public DoubleRefFromArrayElement(final double[] arr, final int off) {
            this.arr = arr;
            this.off = off;
        }
        public DoubleRefFromArrayElement(final $double x, final int off) throws NRException {
            if (!(x instanceof DoubleRefFromArrayElement))
                throw new NRException("!(x instanceof DoubleRefFromArrayElement)");
            DoubleRefFromArrayElement temp = (DoubleRefFromArrayElement) x;
            this.arr = temp.arr;
            this.off = temp.off + off;
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
            StringBuilder sb = new StringBuilder();
            sb.append("{ ");
            for (double x : arr)
                sb.append(x + ", ");
            sb.append("} at " + off);
            return sb.toString();
        }
        @Override
        public double $_(int n) throws NRException {
            if (off + n >= arr.length || off + n < 0)
                throw new NRException("out of range offset used for offset $double");
            return arr[off+n];
        }
        @Override
        public void $_(int n, double v) throws NRException {
            if (off + n >= arr.length || off + n < 0)
                throw new NRException("out of range offset used for offset $double");
            arr[off+n] = v;
        }
	}
	
	static class DoubleRefFrom2DArrayElement implements $double {
        private final double[][] arr;
        private int i, j;
        public DoubleRefFrom2DArrayElement(final double[][] arr, final int i, final int j) {
            this.arr = arr;
            this.i = i;
            this.j = j;
        }
        @Override
        public void $(double t) {
            arr[i][j] = t;
        }
        @Override
        public double $() {
            return arr[i][j];
        }
        @Override
        public String toString() {
            return String.valueOf(arr[i][j]);
        }
        @Override
        public void $_(int n, double v) throws NRException {
            throw new NRException();
            
        }
        @Override
        public double $_(int n) throws NRException {
            throw new NRException();
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
        public T $$() throws NRException {
            return t.copyOut();
        }
        @Override
        public void $$(T t) throws NRException {
            this.t.copyIn(t);
        }
        @Override
        public String toString() {
            return t.toString();
        }
    }
}
