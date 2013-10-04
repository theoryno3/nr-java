
package com.snuggy.nr.util;

import java.util.*;

public class LongIndexArrayOfLong {
    private Map<Long,Long> map = new HashMap<Long,Long>();
    
    public void set(long index, long v) {
        if (index < 0)
            throw new ArrayIndexOutOfBoundsException("index is " + index);
        map.put(index, v);
    }
    public long get(long index) {
        if (!map.containsKey(index))
            throw new ArrayIndexOutOfBoundsException("index is " + index);
        return map.get(index);
    }
}
