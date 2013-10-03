
package com.snuggy.nr.chapter04;

import com.snuggy.nr.util.*;

public class NRf3 implements Func_Doub_To_Doub {
    
	private double xsav,ysav;
	private Func_Doub_Doub_Doub_To_Doub func3d;
	
	public void set_func3d(final Func_Doub_Doub_Doub_To_Doub func3d) {
	    this.func3d = func3d;
	}
	
	public void set_ysav(final double ysav) {
	    this.ysav = ysav;
	}
	
	public double xsav() {
	    return xsav;
	}
	
	public void set_xsav(final double xsav) {
	    this.xsav = xsav;
	}
	
	public double eval(final double z)  { 
		// The integrand f.x;y;z/ evaluated at fixed x and y.
		return func3d.eval(xsav,ysav,z);
	}

}
