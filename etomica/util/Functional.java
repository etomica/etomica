package etomica.util;

import etomica.api.IFunction;

/**
 * Interface for the basic features of a functional, which maps a function onto a double.
 */

public interface Functional {
    
    public double f(IFunction x);
 
}