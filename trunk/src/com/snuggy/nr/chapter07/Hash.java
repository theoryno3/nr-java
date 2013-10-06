
package com.snuggy.nr.chapter07;

import java.util.*;

import com.snuggy.nr.refs.*;
import com.snuggy.nr.util.*;

public class Hash<hfnT extends HfnT<keyT,elT>, keyT, elT> extends Hashtable<hfnT, keyT, elT> {

    // Extend the Hashtable class with storage for elements of type elT, and
    // provide methods for storing, retrieving. and erasing elements. key is
    // passed by address in all methods.
    // using Hashtable<keyT,hfnT>::iget;
    // using Hashtable<keyT,hfnT>::iset;
    // using Hashtable<keyT,hfnT>::ierase;
    private List<elT> els;
    private Class<elT> elementClass;

    public Hash(final Class<hfnT> hashFunctionClass, 
                final Class<elT> elementClass,
                final int nh, final int nm) throws InstantiationException, IllegalAccessException {
        // Same constructor syntax as Hashtable.
        super(hashFunctionClass, nh, nm);
        this.elementClass = elementClass;
        els = new ArrayList<elT>();
        for (int i = 0; i < nm; i++)
            els.add(null);
    }

    public void set(final keyT key, final elT el) throws NRException {
        // Store an element el.
        els.set(iset(key), el);
    }

    public int get(final keyT key, final $<elT> el) throws NRException {
        // Retrieve an element into el. Returns 0 if no element is stored
        // under key, or 1 for success.
        int ll = iget(key);
        if (ll < 0)
            return 0;
        el.$(els.get(ll));
        return 1;
    }

    public elT get(final keyT key) throws NRException, InstantiationException, IllegalAccessException {
        // Store or retrieve an element using subscript notation for its key.
        // Returns a reference that can be used as an l-value.
        int ll = iget(key);
        if (ll < 0) {
            ll = iset(key);
            els.set(ll, elementClass.newInstance());
        }
        return els.get(ll);
    }

    public int count(final keyT key) throws NRException {
        // Return the number of elements stored under key, that is, either 0 or
        // 1.
        int ll = iget(key);
        return (ll < 0 ? 0 : 1);
    }

    public int erase(final keyT key) throws NRException {
        // Erase an element. Returns 1 for success, or 0 if no element is
        // stored under key.
        return (ierase(key) < 0 ? 0 : 1);
    }

}