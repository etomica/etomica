/*
 * History
 * Created on Aug 30, 2004 by kofke
 */
package etomica;

/**
 * Returns iterates from from the childList of a basis atom.  Behavior is set
 * via iterator directive: if no atom is specified there, all pairs formed from the childList
 * are given; otherwise, if an atom is specified there, pairs will be formed from the
 * childList atoms with the basis' child from which the directive atom is descended.  
 */
public final class ApiIntragroup extends AtomsetIteratorAdapter implements
		AtomsetIteratorBasisDependent {

	/**
	 * Constructor makes iterator that must have basis specified and then be 
	 * reset() before iteration.
	 */
	public ApiIntragroup() {
		super(new ApiInnerVariable(
					new AtomIteratorDirectable(),
					new AtomIteratorSequencerList()));
		pairIterator = (ApiInnerVariable)iterator;
		aiOuter = (AtomIteratorDirectable)pairIterator.getOuterIterator();
		aiInner = (AtomIteratorSequencerList)pairIterator.getInnerIterator();
		aiInner.setSkippingFirst(true);
	}

	/* (non-Javadoc)
	 * @see etomica.AtomsetIteratorBasisDependent#setDirective(etomica.IteratorDirective)
	 */
	//TODO consider a setTarget method instead
	public void setTarget(Atom[] targetAtoms) {
		aiOuter.setTarget(targetAtoms);
		doBoth = targetAtoms.atomCount() == 1;
	}
	
	/**
	 * Puts iterator in a state to begin iteration.
	 */
	public void reset() {
		if(!upListNow) {
			aiInner.setIterationDirection(IteratorDirective.UP);
			upListNow = true;
		}
		pairIterator.reset();
		if(doBoth && !pairIterator.hasNext()) {
			aiInner.setIterationDirection(IteratorDirective.DOWN);
			upListNow = false;
			pairIterator.reset();
		}		
	}

	/**
	 * Specifies the parent atom of the iterates. Pairs given by the
	 * iterator will be formed from the childList of the first atom
	 * in the given array. If argument is null or otherwise does not
	 * specify a childList, iterator will not return iterates 
	 * (hasNext false) until a valid basis is specified.  Length of given
	 * array should match the value returned by setBasis, but if it
	 * is greater no error results; only first atom in array is used.
	 */
	public void setBasis(Atom[] atoms) {
		basis = (atoms != null) ? atoms[0] : null;
		aiOuter.setBasis(atoms);
	}
	
	/**
	 * Returns the next atom pair from the iterator.
	 */
	public Atom[] next() {
		Atom[] next = pairIterator.next();
		if(upListNow && doBoth && !pairIterator.hasNext()) {
			aiInner.setIterationDirection(IteratorDirective.DOWN);
			upListNow = false;
			pairIterator.reset();
		}
		return next;
	}
	
	/**
	 * Returns the number of atom pairs the iterator will return if
	 * reset and iterated in its present state.
	 */
	public int size() {
		int n = ((AtomTreeNodeGroup)basis.node).childList.size();
		return doBoth ? n-1 : n*(n-1)/2;
	}
	
	/**
	 * Indicates whether the given atom pair will be among the iterates
	 * given by the iterator if reset in its present state.  True only
	 * if an iterated pair would match the atoms as ordered in the given
	 * array.
	 */
	public boolean contains(Atom[] atoms) {
		if(atoms==null || atoms[0]==null || atoms[1]==null || atoms[0]==atoms[1]) return false;
		AtomList list = ((AtomTreeNodeGroup)basis.node).childList;
		if(doBoth) return (atoms[0] == basis) && list.contains(atoms[1]);
		else return list.contains(atoms[0]) && list.contains(atoms[1]);
	}

	/**
	 * Returns 1, indicating the the array given to the setBasis method
	 * should have only one element.
	 */
	public int basisSize() {
		return 1;
	}
	
	private final ApiInnerVariable pairIterator;//local, specifically typed copy
	private final AtomIteratorDirectable aiOuter;//local, specifically typed copy
	private final AtomIteratorSequencerList aiInner;//local, specifically typed copy
	private boolean doBoth;//indicates if inner iteration should proceed UP and then DOWN from atom
	private boolean upListNow;//indicates if inner iteration is currently in the UP direction
	private Atom basis;//atom most recently specified in setBasis; used by size() and contains(Atom[]) methods
}
