
package com.snuggy.nr.refs;

public interface ByValue<T> {
    T copyOut();
    void copyIn(T t);
}
