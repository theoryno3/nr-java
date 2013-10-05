package com.snuggy.nr.chapter16;

import static com.snuggy.nr.refs.Refs.*;
import static com.snuggy.nr.util.Static.*;

import com.snuggy.nr.refs.*;

public abstract class Phylagglom {

    // Abstract base class for constructing an agglomerative phylogenetic tree.
    protected int n, root, fsroot; // No. of data points, root node, forced
                                   // root.
    private double seqmax, depmax; // Max. values of seq, dep over the tree.
    protected Phylagglomnode[] t; // The tree.

    protected abstract void premin(final double[][] d, final int[] nextp);

    // Function called before minimum search.
    protected abstract double dminfn(final double[][] d, final int i, final int j);

    // Distance function to be minimized
    protected abstract double dbranchfn(final double[][] d, final int i, final int j);

    // Branch length, node i to mother (j is sister).
    protected abstract double dnewfn(final double[][] d, final int k, final int i, final int j, final int ni,
            final int nj);

    // Distance function for newly constructed nodes.
    protected abstract void drootbranchfn(final double[][] d, final int i, final int j, final int ni, final int nj,
            final $double bi, final $double bj);

    // Sets branch lengths to the final root node.
    // Int comancestor(Int leafa, Int leafb); See text discussion of NJ.

    public Phylagglom(final double[][] dist) throws InstantiationException, IllegalAccessException {
        this(dist, -1);
    }

    public Phylagglom(final double[][] dist, final int fsr) throws InstantiationException, IllegalAccessException {
        n = (nrows(dist));
        fsroot = (fsr);
        t = obj_vec(Phylagglomnode.class, 2 * n - 1);
    }

    public void makethetree(final double[][] dist) {
        // Routine that actually constructs the tree, called by the
        // constructor of a derived class.
        int i, j, k, imin = 0, jmin = 0, node, ntask;
        @SuppressWarnings("unused")
        int ncurr;
        double dd, dmin;
        final double[][] d = doub_mat(dist); // Matrix d is initialized with dist.
        final int[] tp = int_vec(n), nextp = int_vec(n), prevp = int_vec(n), tasklist = int_vec(2 * n + 1);
        final double[] tmp = doub_vec(n);
        for (i = 0; i < n; i++) { // Initializations on leaf elements.
            nextp[i] = i + 1; // nextp and prevp are for looping on the distance
            prevp[i] = i - 1; // matrix even as it becomes sparse.
            tp[i] = i; // tp points from a distance matrix row to a tree
            t[i].ldau = t[i].rdau = -1; // element.
            t[i].nel = 1;
        }
        prevp[0] = nextp[n - 1] = -1; // Signifying end of loop.
        ncurr = n;
        for (node = n; node < 2 * n - 2; node++) { // Main loop!
            premin(d, nextp); // Any calculations needed before min finding.
            dmin = 9.99e99;
            for (i = 0; i >= 0; i = nextp[i]) { // Find i; j pair with min
                                                // distance.
                if (tp[i] == fsroot)
                    continue;
                for (j = nextp[i]; j >= 0; j = nextp[j]) {
                    if (tp[j] == fsroot)
                        continue;
                    if ((dd = dminfn(d, i, j)) < dmin) {
                        dmin = dd;
                        imin = i;
                        jmin = j;
                    }
                }
            }
            i = imin;
            j = jmin;
            t[tp[i]].mo = t[tp[j]].mo = node; // Now set properties of the
                                              // parent
            $(t[tp[i]].modist, dbranchfn(d, i, j)); // and children.
            $(t[tp[j]].modist, dbranchfn(d, j, i));
            t[node].ldau = tp[i];
            t[node].rdau = tp[j];
            t[node].nel = t[tp[i]].nel + t[tp[j]].nel;
            for (k = 0; k >= 0; k = nextp[k]) { // Get new-node distances.
                tmp[k] = dnewfn(d, k, i, j, t[tp[i]].nel, t[tp[j]].nel);
            }
            for (k = 0; k >= 0; k = nextp[k])
                d[i][k] = d[k][i] = tmp[k];
            tp[i] = node; // New node replaces child i in dist. matrix, while
                          // child
            if (prevp[j] >= 0)
                nextp[prevp[j]] = nextp[j]; // j gets patched around.
            if (nextp[j] >= 0)
                prevp[nextp[j]] = prevp[j];
            ncurr--;
        } // End of main loop.
        i = 0;
        j = nextp[0]; // Set properties of the root node.
        root = node;
        t[tp[i]].mo = t[tp[j]].mo = t[root].mo = root;
        drootbranchfn(d, i, j, t[tp[i]].nel, t[tp[j]].nel, t[tp[i]].modist, t[tp[j]].modist);
        t[root].ldau = tp[i];
        t[root].rdau = tp[j];
        $(t[root].modist, t[root].dep = 0.);
        t[root].nel = t[tp[i]].nel + t[tp[j]].nel;
        // We now traverse the tree computing seq and dep, hints for where
        // to plot nodes in a two-dimensional representation. See text.
        ntask = 0;
        seqmax = depmax = 0.;
        tasklist[ntask++] = root;
        while (ntask > 0) {
            i = tasklist[--ntask];
            if (i >= 0) { // Meaning, process going down the tree.
                t[i].dep = t[t[i].mo].dep + t[i].modist.$();
                if (t[i].dep > depmax)
                    depmax = t[i].dep;
                if (t[i].ldau < 0) { // A leaf node.
                    t[i].seq = seqmax++;
                } else { // Not a leaf node.
                    tasklist[ntask++] = -i - 1;
                    tasklist[ntask++] = t[i].ldau;
                    tasklist[ntask++] = t[i].rdau;
                }
            } else { // Meaning, process coming up the tree.
                i = -i - 1;
                t[i].seq = 0.5 * (t[t[i].ldau].seq + t[t[i].rdau].seq);
            }
        }
    }

    public int comancestor(final int leafa, final int leafb) {
        // Given the node numbers of two leaves, return the node number of
        // their first common ancestor.
        int i, j;
        for (i = leafa; i != root; i = t[i].mo) {
            for (j = leafb; j != root; j = t[j].mo)
                if (i == j)
                    break;
            if (i == j)
                break;
        }
        return i;
    }

}
