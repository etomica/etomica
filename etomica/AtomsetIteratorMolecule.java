/*
 * History
 * Created on Dec 30, 2004 by kofke
 */
package etomica;


/**
 * Implementations of this interface are used by the potential master as the
 * iterators associated with a potential.  These iterators generate molecules
 * (children of species agents) as their iterates, and must necessarily be
 * phase dependent, directable, and targetable.  This interface defines no new
 * methods, but collects the appropriate interfaces into a single interface.
 */
public interface AtomsetIteratorMolecule extends AtomsetIteratorPhaseDependent,
        AtomsetIteratorDirectable, AtomsetIteratorTargetable {

}
