/*
 * History
 * Created on Aug 12, 2004 by kofke
 */
package etomica;

/**
 * Interface indicating that an atom iterator can determine appropriate
 * atoms for iteration given an arbitrary phase.  
 */
public interface AtomIteratorPhaseDependent extends AtomIterator {

	public void setPhase(Phase phase);
	
}
