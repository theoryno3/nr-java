
package com.snuggy.nr.chapter07;

import java.nio.*;

import com.snuggy.nr.util.*;

public class Hashfn1 {

    // Example of an object encapsulating a hash function for use by the
    // class Hashmap.
    private Ranhash hasher; // The actual hash function.
    private int n; // Size of key in bytes.

    public Hashfn1(int nn) {
        n = (nn);
    } // Constructor just saves key size.

    public long fn(final byte[] key_arr, final int key_off) throws NRException {
        ByteBuffer key = ByteBuffer.wrap(key_arr);
        // Function that returns hash from key.
        // Uint *k;
        // Ullong *kk;
        switch (n) {
        case 4:
            // k = (Uint *)key;
            // return hasher.int64(*k); // Return 64-bit hash of 32-bit key.
            if (key_off % 4 != 0)
                throw new NRException("Want key offset that is multiple of 4.");
            return hasher.int64(key.asIntBuffer().get(key_off / 4));
        case 8:
            // kk = (Ullong *)key;
            // return hasher.int64(*kk); // Return 64-bit hash of 64-bit key.
            if (key_off % 8 != 0)
                throw new NRException("Want key offset that is multiple of 8.");
            return hasher.int64(key.asLongBuffer().get(key_off / 8));
        default:
            throw new NRException("Hashfn1 is for 4 or 8 byte keys only.");
        }
    }

}
