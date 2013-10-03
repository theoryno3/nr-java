
package com.snuggy.nr.chapter07;

public class Hashfn2Size4<keyT,elT> extends Hashfn2<keyT,elT> {
    
    private Hashfn2Size4<keyT,elT> f = new Hashfn2Size4<>();

    protected Hashfn2Size4() {
        super(4);
    }

    @Override
    public long fn(keyT key) {
        return f.fn(key);
    }

    @Override
    public byte[] keyToBytes(keyT key) {
        // TODO Auto-generated method stub
        return null;
    }

}
