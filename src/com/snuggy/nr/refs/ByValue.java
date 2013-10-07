
package com.snuggy.nr.refs;

import com.snuggy.nr.util.*;

public interface ByValue<T> {
    T copyOut() throws NRException;
    void copyIn(T t) throws NRException;
}
