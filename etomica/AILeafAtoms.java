/*
 * History
 * Created on Aug 13, 2004 by kofke
 */
package etomica;

/**
 * @author kofke
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class AILeafAtoms extends AtomIteratorListAdapter {

	/**
	 * 
	 */
	public AILeafAtoms() {
		super(new AtomIteratorListSimple());
		
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
