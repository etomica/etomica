/*
 * History
 * Created on Aug 13, 2004 by kofke
 */
package etomica;

/**
 * @author kofke
 * Iterator that will loop over all leaf atoms in a phase.
 */
public class AILeafAtoms extends AtomIteratorAdapter {

	private final AtomIteratorListSimple listIterator;
	
	/**
	 * Creates iterator with no phase specified. Iteration will
	 * return no atoms until a call to setPhase is performed.
	 */
	public AILeafAtoms() {
		super(new AtomIteratorListSimple());
		listIterator = (AtomIteratorListSimple)iterator;
	}

	public AILeafAtoms(Phase phase){
		this();
		setPhase(phase);
	}
	
	/* (non-Javadoc)
	 * @see etomica.AtomIteratorListAdapter#setPhase()
	 */
	public void setPhase(Phase phase) {
		listIterator.setList(phase.speciesMaster.atomList);
	}

}
