/*
 * History
 * Created on Aug 13, 2004 by kofke
 */
package etomica;

/**
 * @author kofke
 * Iterator that will loop over all leaf atoms in a phase. Can
 * be configured to iterate all leaf atoms, or only those of
 * a particular species.
 */
public final class AtomIteratorLeafAtoms extends AtomIteratorAdapter implements AtomsetIteratorSpeciesDependent {

	/**
	 * Creates iterator with no phase specified. Iteration will
	 * return no atoms until a call to setPhase is performed.
	 */
	public AtomIteratorLeafAtoms() {
		super(new AtomIteratorList());
		listIterator = (AtomIteratorList)iterator;
	}

	/**
	 * Creates iterator conditioned to give all leaf atoms of the
	 * specified phase.  Call to reset() is required before
	 * beginning iteration.
	 */
	public AtomIteratorLeafAtoms(Phase phase){
		this();
		setPhase(phase);
	}
	
	/**
	 * Configures iterator to form its iterates from the
	 * leaf atoms of the given phase.  If a species was
	 * previously (or subsequently) set, iterates will
	 * be the leaf atoms of under the species in the specified
	 * phase.  If null is given, iterator will return no
	 * iterates until a phase is specified via another call
	 * to the method.
	 */
	public void setPhase(Phase phase) {
		this.phase = phase;
		if (phase==null) {
			listIterator.setList(null);
		}
		else {
			listIterator.setList(phase.speciesMaster.atomList);
			if (species != null) {
				listIterator.setFirst(phase.getAgent(species).firstLeafAtomTab);
				listIterator.setTerminator(null);
			}
		}
		unset();
	}
	
	/**
	 * Configures iterator to form its iterates from the leaf
	 * atoms under the specified species in the previously
	 * (or subsequently) set phase.  If null or zero-length
	 * argument is given, no species is set and subsequent
	 * iteration will be over all leaf atoms in phase.
	 */
	public void setSpecies(Species[] species) {
		this.species = (species != null && species.length > 0) ? species[0] : null;
		setPhase(phase);
	}

	private final AtomIteratorList listIterator;
	private Phase phase;
	private Species species;
	
}
