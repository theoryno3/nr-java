
package com.snuggy.nr.chapter17;

import static com.snuggy.nr.chapter17.Static.*;
import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import com.snuggy.nr.chapter02.*;
import com.snuggy.nr.chapter07.*;
import com.snuggy.nr.util.*;

public class Stochsim {

    private static final int mm = 3; // Set number of reactions.
    private static final int nn = 4; // Set number of species.

    // Object for stochastic simulation of a set of chemical reactions.
    private double[] s; // Vector of species numbers.
    private double[] a; // Vector of rates.
    private double[][] instate, outstate;
    private NRsparseCol[] outchg, depend; // Sparse matrices ij and Gij
    private int[] pr; // Priority list.
    private double t, asum;
    private Ran ran;
    
    public double[] s() {
        return s;
    }

    // typedef Doub(Stochsim::*rateptr)(); Obscure C++ used to create a vector
    // dispatch of function pointers to the rate functions.
    private Func_Void_To_Doub dispatch_arr[];
    private int dispatch_off;
    @SuppressWarnings("unused")
    private int rateptr_off;
    // begin user section
    // Replace this section, using as a template the example (17.7.1) shown
    // here,
    // by the particulars of your reaction network. If you have a large number
    // of reactions, you will want to generate the matrices instate and outstate
    // externally, and pass them as globals (or read them here).
    double k0, k1, k2; // Declare any rate constants needed.

    public double rate0() {
        return k0 * s[0] * s[1];
    } // Your rate functions go here.

    public double rate1() {
        return k1 * s[1] * s[2];
    }

    public double rate2() {
        return k2 * s[2];
    }

    public void describereactions() throws NRException {
        // You provide a function with this name that sets any constants that
        // you have defined and sets the instate and outstate matrices to
        // describe your reactions.
        k0 = 0.01;
        k1 = .1;
        k2 = 1.;
        double indat[] = { // The reactant matrix ij .
        1., 0., 0., 1., 1., 0., 0., 1., 1., 0., 0., 0. };
        instate = doub_mat(nn, mm, indat);
        double outdat[] = { // The state change matrix ij .
        -1., 0., 0., 1., -1., 0., 0., 1., -1., 0., 0., 1. };
        outstate = doub_mat(nn, mm, outdat);
        // dispatch[0] = &Stochsim::rate0; // You must also point the dispatch
        // table
        // entries to the correct rate functions.
        dispatch_arr[0] = new Func_Void_To_Doub() {
            @Override
            public double eval() {
                return rate0();
            }
        };
        // dispatch[1] = &Stochsim::rate1;
        dispatch_arr[1] = new Func_Void_To_Doub() {
            @Override
            public double eval() {
                return rate1();
            }
        };
        // dispatch[2] = &Stochsim::rate2;
        dispatch_arr[2] = new Func_Void_To_Doub() {
            @Override
            public double eval() {
                return rate2();
            }
        };
    }

    // end user section

    public Stochsim(double[] sinit) throws NRException, InstantiationException, IllegalAccessException {
        this(sinit, 1);
    }

    public Stochsim(double[] sinit, final int seed) throws NRException, InstantiationException, IllegalAccessException {
        // Constructor. Input initial species numbers and an optional random
        // seed.
        s = (sinit);
        a = doub_arr(mm, 0.);
        outchg = obj_arr(NRsparseCol.class, mm);
        depend = obj_arr(NRsparseCol.class, mm);
        pr = int_arr(mm);
        t = (0.);
        asum = (0.);
        ran = new Ran(seed);
        // dispatch = (new rateptr[mm]);
        dispatch_arr = (new Func_Void_To_Doub[mm]);
        dispatch_off = 0;
        int i, j, k, d;
        describereactions();
        sparmatfill(outchg, outstate);
        double[][] dep = doub_mat(mm, mm);
        for (i = 0; i < mm; i++)
            for (j = 0; j < mm; j++) { // Logical matrix multiply calculates the
                d = 0; // dependency matrix.
                for (k = 0; k < nn; k++)
                    d = ((d != 0) || ((instate[k][i] != 0.0) && (outstate[k][j] != 0.0))) ? 1 : 0;
                dep[i][j] = d;
            }
        sparmatfill(depend, dep);
        for (i = 0; i < mm; i++) { // Calculate all initial rates.
            pr[i] = i;
            a[i] = dispatch_arr[dispatch_off + i].eval();
            asum += a[i];
        }
    }

    // ~Stochsim() {delete [] dispatch;}

    public double step() {
        // Take a single stochastic step (one reaction) and return the new time.
        int i, n, m, k = 0;
        double tau, atarg, sum, anew;
        if (asum == 0.) {
            t *= 2.;
            return t;
        } // Rare: All reactions have stopped
        // exactly, so double the time until the user notices!
        tau = -log(ran.doub()) / asum;
        atarg = ran.doub() * asum;
        sum = a[pr[0]];
        while (sum < atarg)
            sum += a[pr[++k]]; // Equation (17.7.7).
        m = pr[k];
        if (k > 0)
            SWAP(pr, k, k - 1); // Move reaction up on the priority list.
        if (k == mm - 1)
            asum = sum; // Free update of asum fixes accumulated
        // roundoff.
        n = outchg[m].nvals();
        for (i = 0; i < n; i++) { // Apply state change vector.
            k = outchg[m].row_ind()[i];
            s[k] += outchg[m].val()[i];
        }
        n = depend[m].nvals();
        for (i = 0; i < n; i++) { // Recalculate rates required by depen
            k = depend[m].row_ind()[i]; // dency matrix.
            anew = dispatch_arr[dispatch_off + k].eval();
            asum += (anew - a[k]);
            a[k] = anew;
        }
        if (t * asum < 0.1) // Rare: Rates heading toward zero.
            for (asum = 0., i = 0; i < mm; i++)
                asum += a[i]; // Better recalculate asum.
        return (t += tau);
    }
}
