
package com.snuggy.nr.util;

import java.util.*;

public class LongIndexArrayOfInt {
    private Map<Long,Integer> map = new HashMap<Long,Integer>();
    
    public void set(long index, int v) {
        if (index < 0)
            throw new ArrayIndexOutOfBoundsException("index is " + index);
        map.put(index, v);
    }
    public int get(long index) {
        if (!map.containsKey(index))
            throw new ArrayIndexOutOfBoundsException("index is " + index);
        return map.get(index);
    }
}
