/*
 * History
 * Created on Aug 31, 2004 by kofke
 */
package etomica;

/**
 * Iterator that is set using an array of Species.  Extends
 * AtomsetIteratorPhaseDependent because the species-dependent
 * iterates are necessarily phase-dependent too.
 */
public interface AtomsetIteratorSpeciesDependent extends
		AtomsetIteratorPhaseDependent {

	public void setSpecies(Species[] species);
}
