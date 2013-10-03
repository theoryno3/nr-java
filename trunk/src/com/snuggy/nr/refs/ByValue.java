
package com.snuggy.nr.refs;

public interface ByValue<T> {
    void copyIn(T t);
    T copyOut();
}
