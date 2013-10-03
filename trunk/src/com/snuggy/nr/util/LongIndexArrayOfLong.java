
package com.snuggy.nr.util;

import java.util.*;

public class LongIndexArrayOfLong {
    private Map<Long,Long> map = new HashMap<Long,Long>();
    
    public void set(long index, long v) {
        map.put(index, v);
    }
    public long get(long index) {
        return map.get(index);
    }
}
