/*
 * History
 * Created on Aug 30, 2004 by kofke
 */
package etomica;

/**
 * @author kofke
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class AtomIteratorDirectable extends AtomIteratorAdapter implements
		AtomsetIteratorDirectable {

	/**
	 * @param iterator
	 */
	public AtomIteratorDirectable() {
		super(new AtomIteratorListSimple());
		listIterator = (AtomIteratorListSimple)iterator;
	}

	/* (non-Javadoc)
	 * @see etomica.AtomsetIteratorDirectable#setDirective(etomica.IteratorDirective)
	 */
	public void setDirective(IteratorDirective id) {
		targetAtom = id.atom1();
	}

	/* (non-Javadoc)
	 * @see etomica.AtomsetIteratorDirectable#setBasis(etomica.Atom[])
	 */
	public void setBasis(Atom[] atoms) {
		//if targetAtom == null loop over all children
		listIterator.setList(((AtomTreeNodeGroup)atoms[0].node).childList);
		
		//return child of basis that is or is above targetAtom (if in hierarchy of basis)
		//do no looping if not in hierarchy of basis

	}

	/* (non-Javadoc)
	 * @see etomica.AtomsetIteratorDirectable#basisSize()
	 */
	public int basisSize() {
		return 1;
	}

	private AtomIteratorListSimple listIterator;
	private Atom targetAtom;
}
