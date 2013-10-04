
package com.snuggy.nr.chapter07;

import static com.snuggy.nr.util.Static.*;

import com.snuggy.nr.util.*;

public class Hashtable<hfnT extends HfnT<keyT,eltT>, keyT, eltT> {

    // Instantiate a hash table, with methods for maintaining a one-to-one
    // correspondence between arbitrary keys and unique integers in a
    // specified range.
    private int nhash, nmax, nn, ng;
    private LongIndexArrayOfInt htable, next, garbg;
    private LongIndexArrayOfLong thehash;
    private hfnT hash; // An instance of a hash function object.

    // Hashtable(Int nh, Int nv);
    // Constructor. Arguments are size of hash table and max number of stored
    // elements (keys).
    // Int iget(const keyT &key); Return integer for a previously set key.
    // Int iset(const keyT &key); Return unique integer for a new key.
    // Int ierase(const keyT &key); Erase a key.
    // Int ireserve(); Reserve an integer (with no key).
    // Int irelinquish(Int k); Un-reserve an integer.

    public Hashtable(Class<hfnT> hashFunctionClass, final int nh, final int nv) 
            throws InstantiationException, IllegalAccessException {
        // Constructor. Set nhash, the size of the hash table, and nmax, the
        // maximum number of elements (keys) that can be accommodated. Allocate
        // arrays appropriately.
        // hash = (sizeof(keyT));
        hash = hashFunctionClass.newInstance();
        nhash = (nh);
        nmax = (nv);
        nn = (0);
        ng = (0);
	    htable = new LongIndexArrayOfInt();
	    next = new LongIndexArrayOfInt();
	    garbg = new LongIndexArrayOfInt();
        thehash = new LongIndexArrayOfLong();
        for (long j = 0; j < nh; j++) {
            htable.set(j, -1);
        } // Signifies empty.
    }

    public int iget(final keyT key) throws NRException {
        // Returns integer in 0..nmax-1 corresponding to key, or 1 if no such
        // key was previously stored.
        long j, k;
        long pp = hash.fn(key); // Get 64-bit hash
        System.out.println(pp);
        j = umod(pp, nhash); // and map it into the hash table.
        for (k = htable.get(j); k != -1; k = next.get(k)) { // Traverse linked list
                                                    // until an ex
            if (thehash.get(k) == pp) { // act match is found.
                return (int) k;
            }
        }
        return -1; // Key was not previously stored.
    }

    public int iset(final keyT key) throws NRException {
        // Returns integer in 0..nmax-1 that will henceforth correspond to key.
        // If key was previously set, return the same integer as before.
        int j, k;
        long kprev = 0;
        long pp = hash.fn(key); // Get 64-bit hash
        j = (int) umod(pp, nhash); // and map it into the hash table.
        if (htable.get(j) == -1) { // Key not in table. Find a free integer, either
            k = (ng != 0) ? garbg.get(--ng) : nn++; // new or previously erased.
            htable.set(j, k);
        } else { // Key might be in table. Traverse list.
            for (k = htable.get(j); k != -1; k = next.get(k)) {
                if (thehash.get(k) == pp) {
                    return k; // Yes. Return previous value.
                }
                kprev = k;
            }
            k = (ng != 0) ? garbg.get(--ng) : nn++; // No. Get new integer.
            next.set(kprev, k);
        }
        if (k >= nmax)
            throw new NRException("storing too many values");
        thehash.set(k, pp); // Store the key at the new or previous integer.
        next.set(k, -1);
        return k;
    }

    public int ierase(final keyT key) throws NRException {
        // Erase a key, returning the integer in 0..nmax-1 erased, or 1 if
        // the key was not previously set.
        int j, k, kprev;
        long pp = hash.fn(key);
        j = (int) umod(pp, nhash);
        if (htable.get(j) == -1)
            return -1; // Key not previously set.
        kprev = -1;
        for (k = htable.get(j); k != -1; k = next.get(k)) {
            if (thehash.get(k) == pp) { // Found key. Splice linked list around it.
                if (kprev == -1)
                    htable.set(j, next.get(k));
                else
                    next.set(kprev, next.get(k));
                garbg.set(ng++, k); // Add k to garbage stack as an available
                                 // integer.
                return k;
            }
            kprev = k;
        }
        return -1; // Key not previously set.
    }

    public int ireserve() throws NRException {
        // Reserve an integer in 0..nmax-1 so that it will not be used by set(),
        // and return its value.
        int k = (ng != 0) ? garbg.get(--ng) : nn++;
        if (k >= nmax)
            throw new NRException("reserving too many values");
        next.set(k, -2);
        return k;
    }

    public int irelinquish(final int k) {
        // Return to the pool an index previously reserved by reserve(), and
        // return it, or return 1 if it was not previously reserved.
        if (next.get(k) != -2) {
            return -1;
        }
        garbg.set(ng++, k);
        return k;
    }

}
