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
		AtomsetIteratorBasisDependent {

	/**
	 * @param iterator
	 */
	public AtomIteratorDirectable() {
		super(new AtomIteratorListSimple());
		listIterator = (AtomIteratorListSimple)iterator;
	}

	/* (non-Javadoc)
	 * @see etomica.AtomsetIteratorBasisDependent#setDirective(etomica.IteratorDirective)
	 */
	public void setTarget(Atom[] targetAtoms) {
		targetAtom = targetAtoms[0];
	}

	/**
	 * Sets the basis for iteration, such that the childList atoms of 
	 * the first atom in the given array will be subject to iteration.
	 * If argument is null or otherwise does not specify a non-leaf atom, 
	 * iterator will be set to give no iterates until a new basis is specified.
	 */
	public void setBasis(Atom[] atoms) {
//		if(atoms == null || atoms.length < 1 || atoms[0] == null || atoms[0].node.isLeaf()) {
//			listIterator.setList(null);
//			return;
//		}
		try {
		//if targetAtom == null loop over all children
		listIterator.setList(((AtomTreeNodeGroup)atoms[0].node).childList);
		
		//return child of basis that is or is above targetAtom (if in hierarchy of basis)
		//do no looping if not in hierarchy of basis
		} catch(Exception e) {listIterator.setList(null);}

	}

	/* (non-Javadoc)
	 * @see etomica.AtomsetIteratorBasisDependent#basisSize()
	 */
	public int basisSize() {
		return 1;
	}

	private AtomIteratorListSimple listIterator;
	private Atom targetAtom;
}
