
package com.snuggy.nr.chapter15;

import com.snuggy.nr.refs.*;

public class State implements ByValue<State> {
    
    private static long set_count = 0;

    // Worked MCMC example: Structure containing the components of x.
    private double lam1, lam2; // 1 and 2
    private double tc; // tc
    private int k1, k2; // k1 and k2
    private double plog; // Set to logP by Plog, below.

    @Override
    public State copyOut() {
        State r = new State();
	    r.lam1 = lam1;
	    r.lam2 = lam2; 
	    r.tc = tc; 
	    r.k1 = k1;
	    r.k2 = k2; 
	    r.plog = plog; 
	    return r;
    }
    
    @Override
    public void copyIn(State t) {
	    lam1 = t.lam1;
	    lam2 = t.lam2; 
	    tc = t.tc; 
	    k1 = t.k1;
	    k2 = t.k2; 
	    plog = t.plog; 
    }

    public State(final double la1, final double la2, final double t, final int kk1, final int kk2) {
        lam1 = (la1);
        lam2 = (la2);
        tc = (t);
        k1 = (kk1);
        k2 = (kk2);
    }
    
    public State() {
    }

    public double tc() {
        return tc;
    }
    
    public void set_tc(double t) {
        tc = t;
    }
    
    public double lam1() {
        return lam1;
    }

    public void set_lam1(double lam) {
        set_count++;
        if (set_count == 1006)
	        System.out.println("count is " + set_count + ", setting to " + lam);
        lam1 = lam;
    }

    public double lam2() {
        return lam2;
    }

    public void set_lam2(double lam) {
        lam2 = lam;
    }

    public int k1() {
        return k1;
    }

    public void set_k1(int k) {
        k1 = k;
    }

    public int k2() {
        return k2;
    }
    
    public void set_k2(int k) {
        k2 = k;
    }
    
    public double plog() {
        return plog;
    }
    
    public void set_plog(double p) {
        plog = p;
    }
}
