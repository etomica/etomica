/*
 * Created on Jul 23, 2004
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica;


/**
 * Interface for a meter that can return a value given an arbitrary atom.
 */
public interface MeterAtomic {
    public double currentValue(Atom a);
}