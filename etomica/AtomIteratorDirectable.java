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
public final class AtomIteratorDirectable extends AtomIteratorAdapter implements
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
		targetAtom = targetAtoms == null ? null : targetAtoms[0];
		setupIterator();
	}

	/**
	 * Sets the basis for iteration, such that the childList atoms of 
	 * the first atom in the given array will be subject to iteration.
	 * If argument is null or otherwise does not specify a non-leaf atom, 
	 * iterator will be set to give no iterates until a new basis is specified.
	 */
	public void setBasis(Atom[] atoms) {
		if(atoms == null || atoms.length < 1) {// || atoms[0] == null || atoms[0].node.isLeaf()) {
			basis = null;
			listIterator.setList(null);
		} else {
			basis = atoms[0];
			setupIterator();
		}
	}
	
	private void setupIterator() {
		try {
			if(targetAtom == null) {
				listIterator.setList(((AtomTreeNodeGroup)basis.node).childList);
			} else {
				//return child of basis that is or is above targetAtom (if in hierarchy of basis)
				//do no looping if not in hierarchy of basis				
				Atom targetNode = targetAtom.node.childWhereDescendedFrom(basis.node).atom;
				littleList.clear();
				littleList.add(targetNode);
				listIterator.setList(littleList);
			}		
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
	private AtomList littleList = new AtomList();
	private Atom basis;
}
