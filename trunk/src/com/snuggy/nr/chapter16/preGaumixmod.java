
package com.snuggy.nr.chapter16;

import static com.snuggy.nr.util.Static.*;

public abstract class preGaumixmod {
    
	// For nonwizards, this is basically a typedev of Mat_mm as an mmmm matrix. 
    // For wizards, what is going on is that we need to set a static variable 
    // mmstat before defining Mat_mm, and this must happen before the 
    // Gaumixmod constructor is invoked.
	protected static int mmstat = -1;
	
	// struct Mat_mm : MatDoub {Mat_mm() : MatDoub(mmstat,mmstat) {} };
	
	protected preGaumixmod(final int mm) {mmstat = mm;}
	
	protected final double[][] Mat_mm() {
	    return doub_mat(mmstat, mmstat);
	}
}
