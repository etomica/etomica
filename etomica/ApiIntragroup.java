/*
 * History
 * Created on Aug 30, 2004 by kofke
 */
package etomica;

/**
 * Returns iterates from from the childList of a single basis atom.  Behavior is set
 * via iterator directive: if no atom is specified there, all pairs formed from the childList
 * are given; otherwise, if an atom is specified there, pairs will be formed from the
 * childList atoms with the basis' child from which the target atom is descended.  
 */
public final class ApiIntragroup extends AtomsetIteratorAdapter implements
		AtomsetIteratorBasisDependent, AtomsetIteratorDirectable, ApiComposite {

	/**
	 * Constructor makes iterator that must have basis specified and then be 
	 * reset() before iteration.
	 */
	public ApiIntragroup() {
		this(new ApiInnerVariable(
					new AtomIteratorBasis(),
					new AtomIteratorSequencerList()));
        ((AtomIteratorSequencerList)aiInner).setNumToSkip(1);
    }
    
    public ApiIntragroup(ApiComposite pairIterator) {
        super(pairIterator);
		aiOuter = (AtomIteratorBasis)pairIterator.getOuterIterator();
		aiInner = (AtomIteratorSequencerList)pairIterator.getInnerIterator();
	}

	public void setTarget(Atom[] targetAtoms) {
		aiOuter.setTarget(targetAtoms);
        oneTarget = targetAtoms[0] != null && (targetAtoms.length == 1 || targetAtoms[1] == null);
	}
	
	/**
	 * Puts iterator in a state to begin iteration.
	 */
	public void reset() {
        //iterate in the prescribed direction from the target
        if(oneTarget) aiInner.setDirection(direction);
        
        //no target given -- iterate over all pairs
        else aiInner.setDirection(IteratorDirective.UP);
        iterator.reset();
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
     * Specifies the direction, which applies only if iterating pairs
     * with a target atom; otherwise, if all pairs from group are indicated,
     * direction is ignored.
     */
	public void setDirection(IteratorDirective.Direction direction) {
        this.direction = direction;
	}
	
	/**
	 * Returns the number of atom pairs the iterator will return if
	 * reset and iterated in its present state.
	 */
	public int size() {
		int n = ((AtomTreeNodeGroup)basis.node).childList.size();
		return oneTarget ? n-1 : n*(n-1)/2;
	}
	
	/**
	 * Returns 1, indicating the the array given to the setBasis method
	 * should have only one element.
	 */
	public int basisSize() {
		return 1;
	}
	
    public AtomIterator getInnerIterator() {
        return (AtomIterator)aiInner;
    }
    
    public AtomIterator getOuterIterator() {
        return aiOuter;
    }
    
	private final AtomIteratorBasis aiOuter;//local, specifically typed copy
	private final AtomIteratorSequencerList aiInner;//local, specifically typed copy
	private boolean oneTarget;
    private IteratorDirective.Direction direction;
	private Atom basis;//atom most recently specified in setBasis; used by size() and contains(Atom[]) methods
}
