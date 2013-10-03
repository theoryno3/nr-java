
package com.snuggy.nr.refs;

import static com.snuggy.nr.refs.Refs.*;

public class Test2 {

    public static void main(String[] argv) {
        $int i = $(17);
        System.out.println("i before is " + i);
        overwrite(i);
        System.out.println("i after is " + i);
    }
    
    public static void overwrite($int t) {
        t.$(23);
    }
}
