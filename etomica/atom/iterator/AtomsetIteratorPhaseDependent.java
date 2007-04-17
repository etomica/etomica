package etomica.atom.iterator;

import etomica.phase.Phase;

/**
 * Interface indicating that an iterator can determine appropriate
 * atoms for iteration given an arbitrary phase.  
 */
public interface AtomsetIteratorPhaseDependent extends AtomsetIterator {

    /**
     * Sets the Phase to pull iterates from
     * @throws a NullPointerException if the Phase is null
     */
	public void setPhase(Phase phase);
	
}
