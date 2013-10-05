
package com.snuggy.nr.chapter16;

import static com.snuggy.nr.util.Static.*;
import static java.lang.Math.*;

import java.io.*;

import com.snuggy.nr.chapter07.*;
import com.snuggy.nr.util.*;

public class Static {

    public static void markovgen(final double[][] atrans, final int[] out) throws NRException {
        markovgen(atrans, out, 0, 1);
    }

    public static void markovgen(final double[][] atrans, final int[] out, final int istart) throws NRException {
        markovgen(atrans, out, istart, 1);
    }

    public static void markovgen(final double[][] atrans, final int[] out, final int istart, final int seed)
            throws NRException {
        // Generate a realization of an M-state Markov model, given its MM
        // transition matrix atrans. The vector out is filled with integers in
        // the range 0:::M 1. The starting state is the optional argument
        // istart (defaults to 0). seed is an optional argument that sets the
        // seed of the random number generator.
        int i, ilo, ihi, ii, j, m = nrows(atrans), n = out.length;
        final double[][] cum = doub_mat(atrans); // Temporary matrix to hold cumulative
                                       // probabilities.
        double r;
        Ran ran = new Ran(seed); // Use the random number generator Ran.
        if (m != ncols(atrans))
            throw new NRException("transition matrix must be square");
        for (i = 0; i < m; i++) { // Fill cum and die if clearly not a
                                  // transition matrix.
            for (j = 1; j < m; j++)
                cum[i][j] += cum[i][j - 1];
            if (abs(cum[i][m - 1] - 1.) > 0.01)
                throw new NRException("transition matrix rows must sum to 1");
        }
        j = istart; // The current state is kept in j.
        out[0] = j;
        for (ii = 1; ii < n; ii++) { // Main loop.
            r = ran.doub() / cum[j][m - 1]; // Slightly-off normalization gets
                                            // corrected here.
            ilo = 0;
            ihi = m;
            while (ihi - ilo > 1) { // Use bisection to find location among the
                                    // cumu
                i = (ihi + ilo) >> 1; // lative probabilities.
                if (r > cum[j][i - 1])
                    ilo = i;
                else
                    ihi = i;
            }
            out[ii] = j = ilo; // Set new current state.
        }
    }

    public static void newick(final Phylagglom p, final char[][] str, String filename) throws FileNotFoundException {
        // Output a phylogenetic tree p in the “Newick” or “New Hampshire”
        // standard format. Text labels for the leaf nodes are input as the
        // rows of str, each a null terminated string. The output file name
        // is specified by filename.
        // FileOutputStream OUT = fopen(filename,"wb");
        PrintStream OUT = new PrintStream(filename);
        int i, s, ntask = 0, n = p.n, root = p.root;
        final int[] tasklist = int_arr(2 * n + 1);
        tasklist[ntask++] = (1 << 16) + root;
        while (ntask-- > 0) { // Depth-first traversal of the tree.
            s = tasklist[ntask] >> 16; // Code indicating context.
            i = tasklist[ntask] & 0xffff; // Node number to be processed.
            System.out.println("s is " + s + ", i is " + i);
            if (s == 1 || s == 2) { // Left or right dau, down.
                tasklist[ntask++] = ((s + 2) << 16) + p.t[i].mo;
                if (p.t[i].ldau >= 0) {
                    printf(OUT, "(");
                    tasklist[ntask++] = (2 << 16) + p.t[i].rdau;
                    tasklist[ntask++] = (1 << 16) + p.t[i].ldau;
                } else {
                    printf(OUT, "%s:%f", new String(str[i], 0, str[i].length-1), p.t[i].modist.$());
                    flush(OUT);
                }
            } else if (s == 3) {
                if (ntask > 0)
                    printf(OUT, ",\n");
            } // Left dau up.
            else if (s == 4) { // Right dau up.
                if (i == root)
                    printf(OUT, ");\n");
                else
                    printf(OUT, "):%f", p.t[i].modist.$());
            }
        }
        close(OUT);
    }
    
    private static void printf(PrintStream out, String format, Object... objs) {
        //System.out.printf(format, objs);
        out.printf(format, objs);
    }
    
    private static void flush(PrintStream out) {
        //System.out.flush();
        out.flush();
    }
    
    public static void close(PrintStream out) {
        out.close();
    }

}
