/*
 * History
 * Created on Aug 31, 2004 by kofke
 */
package etomica;

/**
 * Atom iterator that is set using an array of Species.  Extends
 * AtomIteratorPhaseDependent because the species-dependent
 * iterates are necessarily phase-dependent too.
 */
public interface AtomIteratorSpeciesDependent extends
		AtomIteratorPhaseDependent {

	public void setSpecies(Species[] species);
}
