package com.snuggy.nr.chapter08;

import static com.snuggy.nr.util.Static.*;

import com.snuggy.nr.util.*;

public class Static {

    public static void sort(final double[] arr) throws NRException {
        sort(arr, -1);
    }

    // Here M is the size of subarrays sorted by straight insertion and NSTACK
    // is the required auxiliary storage.
    private static final int sort_M = 7, sort_NSTACK = 64;

    public static void sort(final double[] arr, final int m) throws NRException {
        // Sort an array arr[0..n-1] into ascending numerical order using the
        // Quicksort algorithm. arr is replaced on output by its sorted
        // rearrangement. Normally, the optional argument m should be omitted,
        // but if it is set to a positive value, then only the first m elements
        // of arr are sorted.
        int i, ir, j, k, jstack = -1, l = 0, n = arr.length;
        double a;
        final int[] istack = int_arr(sort_NSTACK);
        if (m > 0)
            n = MIN(m, n); // Use optional argument.
        ir = n - 1;
        for (;;) { // Insertion sort when subarray small enough.
            if (ir - l < sort_M) {
                for (j = l + 1; j <= ir; j++) {
                    a = arr[j];
                    for (i = j - 1; i >= l; i--) {
                        if (arr[i] <= a)
                            break;
                        arr[i + 1] = arr[i];
                    }
                    arr[i + 1] = a;
                }
                if (jstack < 0)
                    break;
                ir = istack[jstack--]; // Pop stack and begin a new round of
                                       // parti
                l = istack[jstack--]; // tioning.
            } else {
                k = (l + ir) >> 1; // Choose median of left, center, and right
                                   // elements
                // as partitioning element a. Also rearrange so that a[l] 
                // a[l+1]  a[ir].
                SWAP(arr, k, l + 1);
                if (arr[l] > arr[ir]) {
                    SWAP(arr, l, ir);
                }
                if (arr[l + 1] > arr[ir]) {
                    SWAP(arr, l + 1, ir);
                }
                if (arr[l] > arr[l + 1]) {
                    SWAP(arr, l, l + 1);
                }
                i = l + 1; // Initialize pointers for partitioning.
                j = ir;
                a = arr[l + 1]; // Partitioning element.
                for (;;) { // Beginning of innermost loop.
                    do
                        i++;
                    while (arr[i] < a); // Scan up to find element > a.
                    do
                        j--;
                    while (arr[j] > a); // Scan down to find element < a.
                    if (j < i)
                        break; // Pointers crossed. Partitioning complete.
                    SWAP(arr, i, j); // Exchange elements.
                } // End of innermost loop.
                arr[l + 1] = arr[j]; // Insert partitioning element.
                arr[j] = a;
                jstack += 2;
                // Push pointers to larger subarray on stack; process smaller
                // subarray immediately.
                if (jstack >= sort_NSTACK)
                    throw new NRException("NSTACK too small in sort.");
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

    // As usual, you can move any other arrays around at the same time as you
    // sort arr. At the risk of being repetitious:

    private static final int sort2_M = 7, sort2_NSTACK = 64;

    public static <U> void sort2(final double[] arr, final U[] brr) throws NRException {
        // Sort an array arr[0..n-1] into ascending order using Quicksort,
        // while making the corresponding rearrangement of the array
        // brr[0..n-1].
        int i, ir, j, k, jstack = -1, l = 0, n = arr.length;
        double a;
        U b;
        final int[] istack = int_arr(sort2_NSTACK);
        ir = n - 1;
        for (;;) { // Insertion sort when subarray small enough.
            if (ir - l < sort2_M) {
                for (j = l + 1; j <= ir; j++) {
                    a = arr[j];
                    b = brr[j];
                    for (i = j - 1; i >= l; i--) {
                        if (arr[i] <= a)
                            break;
                        arr[i + 1] = arr[i];
                        brr[i + 1] = brr[i];
                    }
                    arr[i + 1] = a;
                    brr[i + 1] = b;
                }
                if (jstack < 0)
                    break;
                ir = istack[jstack--]; // Pop stack and begin a new round of
                                       // parti
                l = istack[jstack--]; // tioning.
            } else {
                k = (l + ir) >> 1; // Choose median of left, center, and right
                                   // elements
                // as partitioning element a. Also rearrange so that a[l] 
                // a[l+1]  a[ir].
                SWAP(arr, k, l + 1);
                SWAP(brr, k, l + 1);
                if (arr[l] > arr[ir]) {
                    SWAP(arr, l, ir);
                    SWAP(brr, l, ir);
                }
                if (arr[l + 1] > arr[ir]) {
                    SWAP(arr, l + 1, ir);
                    SWAP(brr, l + 1, ir);
                }
                if (arr[l] > arr[l + 1]) {
                    SWAP(arr, l, l + 1);
                    SWAP(brr, l, l + 1);
                }
                i = l + 1; // Initialize pointers for partitioning.
                j = ir;
                a = arr[l + 1]; // Partitioning element.
                b = brr[l + 1];
                for (;;) { // Beginning of innermost loop.
                    do
                        i++;
                    while (arr[i] < a); // Scan up to find element > a.
                    do
                        j--;
                    while (arr[j] > a); // Scan down to find element < a.
                    if (j < i)
                        break; // Pointers crossed. Partitioning complete.
                    SWAP(arr, i, j); // Exchange elements of both arrays.
                    SWAP(brr, i, j);
                } // End of innermost loop.
                arr[l + 1] = arr[j]; // Insert partitioning element in both
                                     // arrays.
                arr[j] = a;
                brr[l + 1] = brr[j];
                brr[j] = b;
                jstack += 2;
                // Push pointers to larger subarray on stack; process smaller
                // subarray immediately.
                if (jstack >= sort2_NSTACK)
                    throw new NRException("NSTACK too small in sort2.");
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

    public static void sort(final int[] arr) throws NRException {
        sort(arr, -1);
    }

    // Here M is the size of subarrays sorted by straight insertion and NSTACK
    // is the required auxiliary storage.
    private static final int int_sort_M = 7, int_sort_NSTACK = 64;

    public static void sort(final int[] arr, final int m) throws NRException {
        // Sort an array arr[0..n-1] into ascending numerical order using the
        // Quicksort algorithm. arr is replaced on output by its sorted
        // rearrangement. Normally, the optional argument m should be omitted,
        // but if it is set to a positive value, then only the first m elements
        // of arr are sorted.
        int i, ir, j, k, jstack = -1, l = 0, n = arr.length;
        int a;
        final int[] istack = int_arr(int_sort_NSTACK);
        if (m > 0)
            n = MIN(m, n); // Use optional argument.
        ir = n - 1;
        for (;;) { // Insertion sort when subarray small enough.
            if (ir - l < int_sort_M) {
                for (j = l + 1; j <= ir; j++) {
                    a = arr[j];
                    for (i = j - 1; i >= l; i--) {
                        if (arr[i] <= a)
                            break;
                        arr[i + 1] = arr[i];
                    }
                    arr[i + 1] = a;
                }
                if (jstack < 0)
                    break;
                ir = istack[jstack--]; // Pop stack and begin a new round of
                                       // parti
                l = istack[jstack--]; // tioning.
            } else {
                k = (l + ir) >> 1; // Choose median of left, center, and right
                                   // elements
                // as partitioning element a. Also rearrange so that a[l] 
                // a[l+1]  a[ir].
                SWAP(arr, k, l + 1);
                if (arr[l] > arr[ir]) {
                    SWAP(arr, l, ir);
                }
                if (arr[l + 1] > arr[ir]) {
                    SWAP(arr, l + 1, ir);
                }
                if (arr[l] > arr[l + 1]) {
                    SWAP(arr, l, l + 1);
                }
                i = l + 1; // Initialize pointers for partitioning.
                j = ir;
                a = arr[l + 1]; // Partitioning element.
                for (;;) { // Beginning of innermost loop.
                    do
                        i++;
                    while (arr[i] < a); // Scan up to find element > a.
                    do
                        j--;
                    while (arr[j] > a); // Scan down to find element < a.
                    if (j < i)
                        break; // Pointers crossed. Partitioning complete.
                    SWAP(arr, i, j); // Exchange elements.
                } // End of innermost loop.
                arr[l + 1] = arr[j]; // Insert partitioning element.
                arr[j] = a;
                jstack += 2;
                // Push pointers to larger subarray on stack; process smaller
                // subarray immediately.
                if (jstack >= int_sort_NSTACK)
                    throw new NRException("NSTACK too small in sort.");
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

    // As usual, you can move any other arrays around at the same time as you
    // sort arr. At the risk of being repetitious:

    private static final int int_sort2_M = 7, int_sort2_NSTACK = 64;

    public static <U> void sort2(final int[] arr, final U[] brr) throws NRException {
        // Sort an array arr[0..n-1] into ascending order using Quicksort,
        // while making the corresponding rearrangement of the array
        // brr[0..n-1].
        int i, ir, j, k, jstack = -1, l = 0, n = arr.length;
        int a;
        U b;
        final int[] istack = int_arr(int_sort2_NSTACK);
        ir = n - 1;
        for (;;) { // Insertion sort when subarray small enough.
            if (ir - l < int_sort2_M) {
                for (j = l + 1; j <= ir; j++) {
                    a = arr[j];
                    b = brr[j];
                    for (i = j - 1; i >= l; i--) {
                        if (arr[i] <= a)
                            break;
                        arr[i + 1] = arr[i];
                        brr[i + 1] = brr[i];
                    }
                    arr[i + 1] = a;
                    brr[i + 1] = b;
                }
                if (jstack < 0)
                    break;
                ir = istack[jstack--]; // Pop stack and begin a new round of
                                       // parti
                l = istack[jstack--]; // tioning.
            } else {
                k = (l + ir) >> 1; // Choose median of left, center, and right
                                   // elements
                // as partitioning element a. Also rearrange so that a[l] 
                // a[l+1]  a[ir].
                SWAP(arr, k, l + 1);
                SWAP(brr, k, l + 1);
                if (arr[l] > arr[ir]) {
                    SWAP(arr, l, ir);
                    SWAP(brr, l, ir);
                }
                if (arr[l + 1] > arr[ir]) {
                    SWAP(arr, l + 1, ir);
                    SWAP(brr, l + 1, ir);
                }
                if (arr[l] > arr[l + 1]) {
                    SWAP(arr, l, l + 1);
                    SWAP(brr, l, l + 1);
                }
                i = l + 1; // Initialize pointers for partitioning.
                j = ir;
                a = arr[l + 1]; // Partitioning element.
                b = brr[l + 1];
                for (;;) { // Beginning of innermost loop.
                    do
                        i++;
                    while (arr[i] < a); // Scan up to find element > a.
                    do
                        j--;
                    while (arr[j] > a); // Scan down to find element < a.
                    if (j < i)
                        break; // Pointers crossed. Partitioning complete.
                    SWAP(arr, i, j); // Exchange elements of both arrays.
                    SWAP(brr, i, j);
                } // End of innermost loop.
                arr[l + 1] = arr[j]; // Insert partitioning element in both
                                     // arrays.
                arr[j] = a;
                brr[l + 1] = brr[j];
                brr[j] = b;
                jstack += 2;
                // Push pointers to larger subarray on stack; process smaller
                // subarray immediately.
                if (jstack >= int_sort2_NSTACK)
                    throw new NRException("NSTACK too small in sort2.");
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

    public static <T extends Comparable<T>> void sort(final T[] arr) throws NRException {
        sort(arr, -1);
    }

    // Here M is the size of subarrays sorted by straight insertion and NSTACK
    // is the required auxiliary storage.
    private static final int T_sort_M = 7, T_sort_NSTACK = 64;

    public static <T extends Comparable<T>> void sort(final T[] arr, final int m) throws NRException {
        // Sort an array arr[0..n-1] into ascending numerical order using the
        // Quicksort algorithm. arr is replaced on output by its sorted
        // rearrangement. Normally, the optional argument m should be omitted,
        // but if it is set to a positive value, then only the first m elements
        // of arr are sorted.
        int i, ir, j, k, jstack = -1, l = 0, n = arr.length;
        T a;
        final int[] istack = int_arr(T_sort_NSTACK);
        if (m > 0)
            n = MIN(m, n); // Use optional argument.
        ir = n - 1;
        for (;;) { // Insertion sort when subarray small enough.
            if (ir - l < T_sort_M) {
                for (j = l + 1; j <= ir; j++) {
                    a = arr[j];
                    for (i = j - 1; i >= l; i--) {
                        // if (arr[i] <= a)
                        if (arr[i].compareTo(a) <= 0)
                            break;
                        arr[i + 1] = arr[i];
                    }
                    arr[i + 1] = a;
                }
                if (jstack < 0)
                    break;
                ir = istack[jstack--]; // Pop stack and begin a new round of
                                       // parti
                l = istack[jstack--]; // tioning.
            } else {
                k = (l + ir) >> 1; // Choose median of left, center, and right
                                   // elements
                // as partitioning element a. Also rearrange so that a[l] 
                // a[l+1]  a[ir].
                SWAP(arr, k, l + 1);
                // if (arr[l] > arr[ir]) {
                if (arr[l].compareTo(arr[ir]) > 0) {
                    SWAP(arr, l, ir);
                }
                if (arr[l + 1].compareTo(arr[ir]) > 0) {
                    SWAP(arr, l + 1, ir);
                }
                if (arr[l].compareTo(arr[l + 1]) > 0) {
                    SWAP(arr, l, l + 1);
                }
                i = l + 1; // Initialize pointers for partitioning.
                j = ir;
                a = arr[l + 1]; // Partitioning element.
                for (;;) { // Beginning of innermost loop.
                    do
                        i++;
                    while (arr[i].compareTo(a) < 0); // Scan up to find element
                                                     // > a.
                    do
                        j--;
                    while (arr[j].compareTo(a) > 0); // Scan down to find
                                                     // element < a.
                    if (j < i)
                        break; // Pointers crossed. Partitioning complete.
                    SWAP(arr, i, j); // Exchange elements.
                } // End of innermost loop.
                arr[l + 1] = arr[j]; // Insert partitioning element.
                arr[j] = a;
                jstack += 2;
                // Push pointers to larger subarray on stack; process smaller
                // subarray immediately.
                if (jstack >= T_sort_NSTACK)
                    throw new NRException("NSTACK too small in sort.");
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

    // As usual, you can move any other arrays around at the same time as you
    // sort arr. At the risk of being repetitious:

    private static final int T_sort2_M = 7, T_sort2_NSTACK = 64;

    public static <T extends Comparable<T>, U> void sort2(final T[] arr, final U[] brr) throws NRException {
        // Sort an array arr[0..n-1] into ascending order using Quicksort,
        // while making the corresponding rearrangement of the array
        // brr[0..n-1].
        int i, ir, j, k, jstack = -1, l = 0, n = arr.length;
        T a;
        U b;
        final int[] istack = int_arr(T_sort2_NSTACK);
        ir = n - 1;
        for (;;) { // Insertion sort when subarray small enough.
            if (ir - l < T_sort2_M) {
                for (j = l + 1; j <= ir; j++) {
                    a = arr[j];
                    b = brr[j];
                    for (i = j - 1; i >= l; i--) {
                        if (arr[i].compareTo(a) <= 0)
                            break;
                        arr[i + 1] = arr[i];
                        brr[i + 1] = brr[i];
                    }
                    arr[i + 1] = a;
                    brr[i + 1] = b;
                }
                if (jstack < 0)
                    break;
                ir = istack[jstack--]; // Pop stack and begin a new round of
                                       // parti
                l = istack[jstack--]; // tioning.
            } else {
                k = (l + ir) >> 1; // Choose median of left, center, and right
                                   // elements
                // as partitioning element a. Also rearrange so that a[l] 
                // a[l+1]  a[ir].
                SWAP(arr, k, l + 1);
                SWAP(brr, k, l + 1);
                if (arr[l].compareTo(arr[ir]) > 0) {
                    SWAP(arr, l, ir);
                    SWAP(brr, l, ir);
                }
                if (arr[l + 1].compareTo(arr[ir]) > 0) {
                    SWAP(arr, l + 1, ir);
                    SWAP(brr, l + 1, ir);
                }
                if (arr[l].compareTo(arr[l + 1]) > 0) {
                    SWAP(arr, l, l + 1);
                    SWAP(brr, l, l + 1);
                }
                i = l + 1; // Initialize pointers for partitioning.
                j = ir;
                a = arr[l + 1]; // Partitioning element.
                b = brr[l + 1];
                for (;;) { // Beginning of innermost loop.
                    do
                        i++;
                    while (arr[i].compareTo(a) < 0); // Scan up to find element
                                                     // > a.
                    do
                        j--;
                    while (arr[j].compareTo(a) > 0); // Scan down to find
                                                     // element < a.
                    if (j < i)
                        break; // Pointers crossed. Partitioning complete.
                    SWAP(arr, i, j); // Exchange elements of both arrays.
                    SWAP(brr, i, j);
                } // End of innermost loop.
                arr[l + 1] = arr[j]; // Insert partitioning element in both
                                     // arrays.
                arr[j] = a;
                brr[l + 1] = brr[j];
                brr[j] = b;
                jstack += 2;
                // Push pointers to larger subarray on stack; process smaller
                // subarray immediately.
                if (jstack >= T_sort2_NSTACK)
                    throw new NRException("NSTACK too small in sort2.");
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

    private static final int dd_sort2_M = 7, dd_sort2_NSTACK = 64;

    public static void sort2(final double[] arr, final double[] brr) throws NRException {
        // Sort an array arr[0..n-1] into ascending order using Quicksort,
        // while making the corresponding rearrangement of the array
        // brr[0..n-1].
        int i, ir, j, k, jstack = -1, l = 0, n = arr.length;
        double a;
        double b;
        final int[] istack = int_arr(dd_sort2_NSTACK);
        ir = n - 1;
        for (;;) { // Insertion sort when subarray small enough.
            if (ir - l < dd_sort2_M) {
                for (j = l + 1; j <= ir; j++) {
                    a = arr[j];
                    b = brr[j];
                    for (i = j - 1; i >= l; i--) {
                        if (arr[i] <= a)
                            break;
                        arr[i + 1] = arr[i];
                        brr[i + 1] = brr[i];
                    }
                    arr[i + 1] = a;
                    brr[i + 1] = b;
                }
                if (jstack < 0)
                    break;
                ir = istack[jstack--]; // Pop stack and begin a new round of
                                       // parti
                l = istack[jstack--]; // tioning.
            } else {
                k = (l + ir) >> 1; // Choose median of left, center, and right
                                   // elements
                // as partitioning element a. Also rearrange so that a[l] 
                // a[l+1]  a[ir].
                SWAP(arr, k, l + 1);
                SWAP(brr, k, l + 1);
                if (arr[l] > arr[ir]) {
                    SWAP(arr, l, ir);
                    SWAP(brr, l, ir);
                }
                if (arr[l + 1] > arr[ir]) {
                    SWAP(arr, l + 1, ir);
                    SWAP(brr, l + 1, ir);
                }
                if (arr[l] > arr[l + 1]) {
                    SWAP(arr, l, l + 1);
                    SWAP(brr, l, l + 1);
                }
                i = l + 1; // Initialize pointers for partitioning.
                j = ir;
                a = arr[l + 1]; // Partitioning element.
                b = brr[l + 1];
                for (;;) { // Beginning of innermost loop.
                    do
                        i++;
                    while (arr[i] < a); // Scan up to find element > a.
                    do
                        j--;
                    while (arr[j] > a); // Scan down to find element < a.
                    if (j < i)
                        break; // Pointers crossed. Partitioning complete.
                    SWAP(arr, i, j); // Exchange elements of both arrays.
                    SWAP(brr, i, j);
                } // End of innermost loop.
                arr[l + 1] = arr[j]; // Insert partitioning element in both
                                     // arrays.
                arr[j] = a;
                brr[l + 1] = brr[j];
                brr[j] = b;
                jstack += 2;
                // Push pointers to larger subarray on stack; process smaller
                // subarray immediately.
                if (jstack >= dd_sort2_NSTACK)
                    throw new NRException("NSTACK too small in sort2.");
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

    public static double select(final int k, final double[] arr) {
        int i, ir, j, l, mid, n = arr.length;
        double a;
        l = 0;
        ir = n - 1;
        for (;;) {
            if (ir <= l + 1) {
                if (ir == l + 1 && arr[ir] < arr[l])
                    SWAP(arr, l, ir);
                return arr[k];
            } else {
                mid = (l + ir) >> 1;
                SWAP(arr, mid, l + 1);
                if (arr[l] > arr[ir])
                    SWAP(arr, l, ir);
                if (arr[l + 1] > arr[ir])
                    SWAP(arr, l + 1, ir);
                if (arr[l] > arr[l + 1])
                    SWAP(arr, l, l + 1);
                i = l + 1;
                j = ir;
                a = arr[l + 1];
                for (;;) {
                    do
                        i++;
                    while (arr[i] < a);
                    do
                        j--;
                    while (arr[j] > a);
                    if (j < i)
                        break;
                    SWAP(arr, i, j);
                }
                arr[l + 1] = arr[j];
                arr[j] = a;
                if (j >= k)
                    ir = j - 1;
                if (j <= k)
                    l = i;
            }
        }
    }

}
