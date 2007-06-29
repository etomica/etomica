/*
 * History
 * Created on Dec 30, 2004 by kofke
 */
package etomica.atom.iterator;



/**
 * This is an interface for iterators that are box dependent, directable, 
 * and targetable.  This interface defines no new methods, but collects the 
 * appropriate interfaces into a single interface.
 */
public interface AtomsetIteratorPDT extends AtomsetIteratorBoxDependent,
        AtomsetIteratorDirectable, AtomsetIteratorTargetable {

}
