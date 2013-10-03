
package com.snuggy.nr.util;

import java.util.*;

public class LongIndexArrayOfInt {
    private Map<Long,Integer> map = new HashMap<Long,Integer>();
    
    public void set(long index, int v) {
        map.put(index, v);
    }
    public int get(long index) {
        return map.get(index);
    }
}
