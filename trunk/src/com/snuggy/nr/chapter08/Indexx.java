
package com.snuggy.nr.chapter08;

import static com.snuggy.nr.util.Static.*;

import com.snuggy.nr.util.*;

public class Indexx<T extends Comparable<T>> {

    private int n;
    private final int[] indx;

    public Indexx(final T[] arr) throws NRException {
        // Constructor. Calls index and stores an index to the array
        // arr[0..n-1].
        indx = int_arr(arr.length);
        index(arr, 0, arr.length);
    }

    public Indexx() {
        indx = null;
    } // Empty constructor. See text.

    public void sort(Class<T> type, final T[] brr) throws NRException, InstantiationException, IllegalAccessException {
        // Sort an array brr[0..n-1] into the order defined by the stored index.
        // brr is replaced on output by its sorted rearrangement.
        if (brr.length != n)
            throw new NRException("bad size in Index sort");
        T[] tmp = obj_arr(type, brr);
        for (int j = 0; j < n; j++)
            brr[j] = tmp[indx[j]];
    }

    public T el(final T[] brr, final int j) {
        // This function, and the next, return the element of brr that would
        // be in sorted position j according to the stored index. The vector
        // brr is not changed.
        return brr[indx[j]];
    }

    // public T el(final T[] brr, final int j) {
    // Same, but return an l-value.
    // return brr[indx[j]];
    // }

    // public void index(final T[] T *arr_arr, final int arr_off, final int nn);
    // This does the actual work of indexing. Normally not called directly by
    // the user, but see
    // text for exceptions.

    public void rank(final int[] irank_ref[]) {
        // Returns a rank table, whose jth element is the rank of arr[j],
        // where arr is the vector originally indexed. The smallest arr[j]
        // has rank 0.
        irank_ref[0] = int_arr(n);
        for (int j = 0; j < n; j++)
            irank_ref[0][indx[j]] = j;
    }

    public void index(final T[] arr_arr, final int arr_off, final int nn) throws NRException {
        // Indexes an array arr[0..nn-1], i.e., resizes and sets indx[0..nn-1]
        // such that arr[indx[j]] is in ascending order
        // for j D 0;1; : : : ;nn-1. Also sets member value n. The input
        // array arr is not changed.
        final int M = 7, NSTACK = 64;
        int i, indxt, ir, j, k, jstack = -1, l = 0;
        T a;
        final int[] istack = int_arr(NSTACK);
        n = nn;
        //indx = int_arr(n);
        ir = n - 1;
        for (j = 0; j < n; j++)
            indx[j] = j;
        for (;;) {
            if (ir - l < M) {
                for (j = l + 1; j <= ir; j++) {
                    indxt = indx[j];
                    a = arr_arr[arr_off + indxt];
                    for (i = j - 1; i >= l; i--) {
                        if (arr_arr[arr_off + indx[i]].compareTo(a) <= 0)
                            break;
                        indx[i + 1] = indx[i];
                    }
                    indx[i + 1] = indxt;
                }
                if (jstack < 0)
                    break;
                ir = istack[jstack--];
                l = istack[jstack--];
            } else {
                k = (l + ir) >> 1;
                SWAP(indx, k, l + 1);
                if (arr_arr[arr_off + indx[l]].compareTo(arr_arr[arr_off + indx[ir]]) > 0) {
                    SWAP(indx, l, ir);
                }
                if (arr_arr[arr_off + indx[l + 1]].compareTo(arr_arr[arr_off + indx[ir]]) > 0) {
                    SWAP(indx, l + 1, ir);
                }
                if (arr_arr[arr_off + indx[l]].compareTo(arr_arr[arr_off + indx[l + 1]]) > 0) {
                    SWAP(indx, l, l + 1);
                }
                i = l + 1;
                j = ir;
                indxt = indx[l + 1];
                a = arr_arr[arr_off + indxt];
                for (;;) {
                    do
                        i++;
                    while (arr_arr[arr_off + indx[i]].compareTo(a) < 0);
                    do
                        j--;
                    while (arr_arr[arr_off + indx[j]].compareTo(a) > 0);
                    if (j < i)
                        break;
                    SWAP(indx, i, j);
                }
                indx[l + 1] = indx[j];
                indx[j] = indxt;
                jstack += 2;
                if (jstack >= NSTACK)
                    throw new NRException("NSTACK too small in index.");
                if (ir - i + 1 >= j - l) {
                    istack[jstack] = ir;
                    istack[jstack - 1] = i;
                    ir = j - 1;
                } else {
                    istack[jstack] = j - 1;
                    istack[jstack - 1] = l;
                    l = i;
                }
            }
        }
    }

}
