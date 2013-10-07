
package com.snuggy.nr.refs;

import static com.snuggy.nr.refs.Refs.*;

import com.snuggy.nr.util.*;

public class TestRefs {

    public static void main(String[] argv) throws InstantiationException, IllegalAccessException, NRException {
	    $<Foo> ref_foo = $(new Foo("foo1"));
        System.out.println("Name before is " + ref_foo);
        probe(ref_foo);
        System.out.println("Name after is " + ref_foo);
        
        $$<Foo> val_foo = $$(new Foo("foo1"));
        System.out.println("Name before is " + val_foo);
        overwrite(val_foo);
        System.out.println("Name after is " + val_foo);
        
        $int i = $(17);
        System.out.println("i before is " + i);
        overwrite(i);
        System.out.println("i after is " + i);
    }
    
    public static void probe($<Foo> foo) {
        System.out.println("probing " + foo.$());
    }
    
    public static void overwrite($$<Foo> foo) throws NRException {
        System.out.println("overwriting " + foo.$$().toString());
        foo.$$(new Foo("foo2"));
    }
    
    public static void overwrite($int t) {
        System.out.println("overwriting " + t.$());
        t.$(23);
    }
}

class Foo implements ByValue<Foo> {
    private String name;
    public Foo() {
        this.name = "no-arg";
    }
    public Foo(String name) {
        this.name = name;
    }
    @Override
    public String toString() {
        return name;
    }
    @Override
    public Foo copyOut() {
        return new Foo(name);
    }
    @Override
    public void copyIn(Foo foo) {
        name = foo.name;
    }
}