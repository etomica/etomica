/*
 * History
 * Created on Aug 30, 2004 by kofke
 */
package etomica;

/**
 * Basic atom iterator that yields the child atoms of a basis that satisfy
 * a target specification (if given one).
 */
public final class AtomIteratorBasis extends AtomIteratorAdapter implements
		AtomsetIteratorBasisDependent {

	/**
	 * Constructor makes iterator in an unset condition; must call reset
	 * before beginning iteration.
	 */
	public AtomIteratorBasis() {
		super(new AtomIteratorListSimple());
		listIterator = (AtomIteratorListSimple)iterator;
	}

	/**
	 * Method to specify a target atom, causing the iterator to yield
	 * only that atom, or its parent in the hierarchy that is a child 
	 * of the current basis.  Iterator yields no iterates if the specified 
	 * target is not in the hierarchy below the basis at the time of reset.
	 * Specifying a null or zero-length array causes releases any target
	 * restrictions, and specifies that the iterator should give all of the
	 * child atoms of the basis.  Only first atom in given array is relevant.
	 * Call to this method leaves iterator unset; call to reset is required
	 * before beginning iteration.
	 */
	public void setTarget(Atom[] targetAtoms) {
		targetAtom = (targetAtoms == null || targetAtoms.length == 0) ? null : targetAtoms[0];
		needSetupIterator = (basis != null);//flag to setup iterator only if presently has a non-null basis
		listIterator.unset();
	}

	/**
	 * Sets the basis for iteration, such that the childList atoms of 
	 * the first atom in the given array will be subject to iteration (within
	 * any specifications given by a prior or subsequent call to setTarget).
	 * If argument is null or otherwise does not specify a non-leaf atom, 
	 * iterator will be conditioned to give no iterates until a new basis 
	 * is specified.  Any atoms beyond the first one in the given array are ignored.
	 */
	public void setBasis(Atom[] atoms) {
		if(atoms == null || atoms.length == 0) {// || atoms[0] == null || atoms[0].node.isLeaf()) {//eliminate last two conditions and allow null basis
			basis = null;
			listIterator.setList(null);
			needSetupIterator = false;
		} else {
			basis = atoms[0];
			needSetupIterator = true;
			listIterator.unset();
		}
	}
	
	/**
	 * Puts iterator in a state ready to begin iteration.
	 */
	public void reset() {
		if(needSetupIterator) setupIterator();
		listIterator.reset();
	}
	
	/**
	 * Common method to complete tasks needed to adjust to new target or basis.
	 * Any call to setBasis or setTarget sets flag that indicates this method should
	 * be invoked upon reset.
	 */
	private void setupIterator() {
		needSetupIterator = false;
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
		} catch(Exception e) {listIterator.setList(null);}//this could happen if childWhereDescendedFrom returns null

	}

	/**
	 * Return 1, indicating that only a single-atom basis is appropriate.
	 */
	public int basisSize() {
		return 1;
	}

	private final AtomIteratorListSimple listIterator;//the wrapped iterator
	private final AtomList littleList = new AtomList();//used to form a list of one iterate if target is specified
	private Atom targetAtom;
	private Atom basis;
	private boolean needSetupIterator = true;//flag to indicate if setupIterator must be called upon reset
}
