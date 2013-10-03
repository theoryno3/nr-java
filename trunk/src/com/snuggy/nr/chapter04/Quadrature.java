
package com.snuggy.nr.chapter04;

import com.snuggy.nr.util.*;

public abstract class Quadrature {

    // Abstract base class for elementary quadrature algorithms.

    protected int n; // Current level of refinement.

    public abstract double next() throws NRException;

    // Returns the value of the integral at the nth stage of refinement.
    // The function next() must be defined in the derived class.

}
