
package com.snuggy.nr.chapter08;

import static com.snuggy.nr.util.Static.*;
import static com.snuggy.nr.chapter08.Static.*;

import com.snuggy.nr.util.*;

public class Heapselect {

    // Object for tracking the m largest values seen thus far in a stream of
    // values.
    private int m, n, srtd;
    private double[] heap;

    public Heapselect(int mm) {
        m = (mm);
        n = (0);
        srtd = (0);
        heap = doub_vec(mm, 1.e99);
    }

    // Constructor. The argument is the number of largest values to track.

    public void add(double val) throws NRException {
        // Assimilate a new value from the stream.
        int j, k;
        if (n < m) { // Heap not yet filled.
            heap[n++] = val;
            if (n == m)
                sort(heap); // Create initial heap by overkill!
        } else {
            if (val > heap[0]) { // Put it on the heap?
                heap[0] = val;
                for (j = 0;;) { // Sift down.
                    k = (j << 1) + 1;
                    if (k > m - 1)
                        break;
                    if (k != (m - 1) && heap[k] > heap[k + 1])
                        k++;
                    if (heap[j] <= heap[k])
                        break;
                    SWAP(heap, k, j);
                    j = k;
                }
            }
            n++;
        }
        srtd = 0; // Mark heap as “unsorted”.
    }

    public double report(int k) throws NRException {
        // Return the kth largest value seen so far. k=0 returns the largest
        // value seen, k=1 the second largest, : : : , k=m-1 the last position
        // tracked. Also, k must be less than the number of previous values
        // assimilated.
        int mm = MIN(n, m);
        if (k > mm - 1)
            throw new NRException("Heapselect k too big");
        if (k == m - 1)
            return heap[0]; // Always free, since top of heap.
        if (!(srtd != 0)) {
            sort(heap);
            srtd = 1;
        } // Otherwise, need to sort the heap.
        return heap[mm - 1 - k];
    }

}
