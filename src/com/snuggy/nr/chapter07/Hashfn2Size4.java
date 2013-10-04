
package com.snuggy.nr.chapter07;

import java.nio.*;

public class Hashfn2Size4<keyT,elT> extends Hashfn2<keyT,elT> {
    
    protected Hashfn2Size4() {
        super(4);
    }

    @Override
    public byte[] keyToBytes(keyT key) {
        ByteBuffer b = ByteBuffer.allocate(4);
        b.putInt(key.hashCode());
        return b.array();
    }
}
