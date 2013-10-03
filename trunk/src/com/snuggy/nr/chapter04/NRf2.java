
package com.snuggy.nr.chapter04;

import static com.snuggy.nr.chapter04.Static.*;

import com.snuggy.nr.util.*;

public class NRf2 implements Func_Doub_To_Doub {
	private NRf3 f3;
	private Func_Doub_Doub_To_Doub z1;
	private Func_Doub_Doub_To_Doub z2;
	
	public NRf2(final Func_Doub_Doub_To_Doub zz1, final Func_Doub_Doub_To_Doub zz2) {
	    f3 = new NRf3();
	    z1 = (zz1);
	    z2 = (zz2);
    }
	
	public NRf3 f3() {
	    return f3;
	}
	
	public double eval(final double y) throws NRException // This is G of eq. (4.8.4).
	{
		f3.set_ysav(y);
		return qgaus(f3,z1.eval(f3.xsav(),y),z2.eval(f3.xsav(),y));
	}
}
