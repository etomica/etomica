/*
 * History
 * Created on Aug 13, 2004 by kofke
 */
package etomica;

/**
 * @author kofke
 * Iterator that will loop over all leaf atoms in a phase.
 */
public class AtomIteratorLeafAtoms extends AtomIteratorAdapter implements AtomIteratorPhaseDependent {

	private final AtomIteratorList listIterator;
	private Phase phase;
	private Species species;
	
	/**
	 * Creates iterator with no phase specified. Iteration will
	 * return no atoms until a call to setPhase is performed.
	 */
	public AtomIteratorLeafAtoms() {
		super(new AtomIteratorList());
		listIterator = (AtomIteratorList)iterator;
	}

	public AtomIteratorLeafAtoms(Phase phase){
		this();
		setPhase(phase);
	}
	
	/* (non-Javadoc)
	 * @see etomica.AtomIteratorListAdapter#setPhase()
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
	}
	
	public void setSpecies(Species species) {
		this.species = species;
		setPhase(phase);
	}

}
