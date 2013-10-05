package com.snuggy.nr.chapter07;


import static com.snuggy.nr.util.Static.*;

import com.snuggy.nr.util.*;

public class Mhash<hfnT extends HfnT<keyT, elT>, keyT, elT> extends Hashtable<hfnT, keyT, elT> {

    // Extend the Hashtable class with storage for elements of type elT,
    // allowing more than one element to be stored under a single key.
    // using Hashtable<keyT,hfnT>::iget;
    // using Hashtable<keyT,hfnT>::iset;
    // using Hashtable<keyT,hfnT>::ierase;
    // using Hashtable<keyT,hfnT>::ireserve;
    // using Hashtable<keyT,hfnT>::irelinquish;
    private elT[] els;
    private final int[] nextsis; // Links to next sister element under a single key.
    private int nextget;

    // Mhash(Int nh, Int nm); Same constructor syntax as Hashtable.
    // Int store(const keyT &key, const elT &el); Store an element under key.
    // Int erase(const keyT &key, const elT &el); Erase a specified element
    // under key.
    // Int count(const keyT &key); Count elements stored under key.
    // Int getinit(const keyT &key); Prepare to retrieve elements from key.
    // Int getnext(elT &el); Retrieve next element specified by getinit.

    public Mhash(Class<hfnT> hashFuncClass, Class<elT> elementClass, final int nh, final int nm) throws InstantiationException, IllegalAccessException {
        super(hashFuncClass, nh, nm);
        nextget = (-1);
        // els(nm);
        els = obj_vec(elementClass, nm);
        nextsis = int_vec(nm);
        for (int j = 0; j < nm; j++) {
            nextsis[j] = -2;
        } // Initialize to “empty”.
    }

    public int store(final keyT key, final elT el) throws NRException {
        // Store an element el under key. Return index in 0..nmax-1, giving
        // the storage location utilized.
        int j, k;
        j = iset(key); // Find root index for this key.
        if (nextsis[j] == -2) { // It is the first object with this key.
            els[j] = el;
            nextsis[j] = -1; // 1 means it is the terminal element.
            return j;
        } else {
            while (nextsis[j] != -1) {
                j = nextsis[j];
            } // Traverse the tree.
            k = ireserve(); // Get a new index and link it into the list.
            els[k] = el;
            nextsis[j] = k;
            nextsis[k] = -1;
            return k;
        }
    }

    public int erase(final keyT key, final elT el) throws NRException {
        // Erase an element el previously stored under key. Return 1 for
        // success, or 0 if no matching element is found. Note: The ==
        // operation must be defined for the type elT.
        int j = -1, kp = -1, kpp = -1;
        int k = iget(key);
        while (k >= 0) {
            if (j < 0 && el.equals(els[k]))
                j = k; // Save index of matching el as j.
            kpp = kp;
            kp = k;
            k = nextsis[k];
        }
        if (j < 0)
            return 0; // No matching el found.
        if (kpp < 0) { // The element el was unique under key.
            ierase(key);
            nextsis[j] = -2;
        } else { // Patch the list.
            if (j != kp)
                els[j] = els[kp]; // Overwrite j with the terminal
                                  // element
            nextsis[kpp] = -1; // and then shorten the list.
            irelinquish(kp);
            nextsis[kp] = -2;
        }
        return 1; // Success.
    }

    public int count(final keyT key) throws NRException {
        // Return the number of elements stored under key, 0 if none.
        int next, n = 1;
        if ((next = iget(key)) < 0)
            return 0;
        while ((next = nextsis[next]) >= 0) {
            n++;
        }
        return n;
    }

    public int getinit(final keyT key) throws NRException {
        // Initialize nextget so that it points to the first element stored
        // under key. Return 1 for success, or 0 if no such element.
        nextget = iget(key);
        return ((nextget < 0) ? 0 : 1);
    }

    public int getnext(final elT el_ref[]) {
        // If nextget points validly, copy its element into el, update nextget
        // to the next element with the same key, and return 1. Otherwise, do
        // not modify el, and return 0.
        if (nextget < 0) {
            return 0;
        }
        el_ref[0] = els[nextget];
        nextget = nextsis[nextget];
        return 1;
    }

}
