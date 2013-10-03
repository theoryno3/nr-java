

package com.snuggy.nr.chapter07;

import static com.snuggy.nr.util.Static.*;

import java.nio.*;

import com.snuggy.nr.util.*;

public abstract class Hashfn2<keyT,elT> implements HfnT<keyT,elT> {

    private static long[] hashfn_tab = long_arr(256); // Defines storage for the
                                                      // lookup table.

    // Another example of an object encapsulating a hash function, allowing
    // arbitrary fixed key sizes or variable-length null terminated strings.
    // The hash function algorithm is self-contained.
    private long h;
    private int n; // Size of key in bytes, when fixed size.

    @SuppressWarnings("unused")
    private Hashfn2() throws NRException {
        throw new NRException("Need a key size");
    }

    protected Hashfn2(final int nn) {
        n = (nn);
        if (n == 1)
            n = 0; // Null terminated string key signaled by n D 0
        h = 0x544B2FBACAAF1684L; // or 1.
        for (int j = 0; j < 256; j++) { // Length 256 lookup table is
                                        // initialized with
            // values from a 64-bit Marsaglia generator
            // stepped 31 times between each.
            for (int i = 0; i < 31; i++) {
                h = (h >>> 7) ^ h;
                h = (h << 11) ^ h;
                h = (h >>> 10) ^ h;
            }
            hashfn_tab[j] = h;
        }
    }
    
    @Override 
    public int keySize() {
        return n;
    }

    @Override 
    public long fn(final keyT key) { // Function that
        // returns hash from key.
        int j;
        // char *k = (char *)key; // Cast the key pointer to char pointer.
        byte[] key_arr_bytes = 
            ByteBuffer.allocate(keySize()).put(keyToBytes(key)).array();
        int[] key_arr = new int[key_arr_bytes.length];
        for (int m = 0; m < key_arr_bytes.length; m++)
            key_arr[m] = ((int) key_arr_bytes[m]) & 0xFF;
        int key_off = 0;
        h = 0xBB40E64DA205B064L;
        j = 0;
        // while ((n != 0) ? (j++ < n) : (*k)) { // Fixed length or else until
        // null.
        while ((n != 0) ? (j++ < n) : (key_arr[key_off] != 0)) { // Fixed
                                                                 // length
                                                                 // or else
                                                                 // until
                                                                 // null.
            h = (h * 7664345821815920749L) ^ hashfn_tab[key_arr[key_off]];
            // k++;
            key_off++;
        }
        return h;
    }

}
