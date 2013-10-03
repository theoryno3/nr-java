
package com.snuggy.nr.chapter07;

public interface HfnT<keyT,elT> {
    long fn(final keyT key);
    int keySize();
    byte[] keyToBytes(final keyT key);
}
