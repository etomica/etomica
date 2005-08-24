/*
 * History
 * Created on Aug 12, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.Phase;
import etomica.atom.AtomsetIterator;

/**
 * Interface indicating that an iterator can determine appropriate
 * atoms for iteration given an arbitrary phase.  
 */
public interface AtomsetIteratorPhaseDependent extends AtomsetIterator {

	public void setPhase(Phase phase);
	
}
