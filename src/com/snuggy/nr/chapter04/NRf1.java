
package com.snuggy.nr.chapter04;

import static com.snuggy.nr.chapter04.Static.*;

import com.snuggy.nr.util.*;

public class NRf1 implements Func_Doub_To_Doub {
    
	private Func_Doub_To_Doub y1;
	private Func_Doub_To_Doub y2;
	private NRf2 f2;
	
	public NRf1(final Func_Doub_To_Doub yy1, final Func_Doub_To_Doub yy2, final Func_Doub_Doub_To_Doub z1, final Func_Doub_Doub_To_Doub z2) {
	    y1 = (yy1);
	    y2 = (yy2);
	    f2 = new NRf2(z1,z2);
    }
	
	public NRf2 f2() {
	    return f2;
	}
	
	public double eval(final double x) throws NRException  {
	// This is H of eq. (4.8.5).
		f2.f3().set_xsav(x);
		return qgaus(f2,y1.eval(x),y2.eval(x));
	}

}
