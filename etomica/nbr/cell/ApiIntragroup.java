/*
 * History
 * Created on Aug 30, 2004 by kofke
 */
package etomica.nbr.cell;

import etomica.ApiInnerVariable;
import etomica.Atom;
import etomica.AtomIteratorBasis;
import etomica.AtomIteratorSequencerList;
import etomica.AtomList;
import etomica.AtomTreeNodeGroup;
import etomica.AtomsetIteratorAdapter;
import etomica.AtomsetIteratorBasisDependent;
import etomica.AtomsetIteratorDirectable;
import etomica.IteratorDirective;
import etomica.Species;

/**
 * Returns iterates from from the childList of a single basis atom.  Behavior is set
 * via iterator directive: if no atom is specified there, all pairs formed from the childList
 * are given; otherwise, if an atom is specified there, pairs will be formed from the
 * childList atoms with the basis' child from which the directive atom is descended.  
 */
public final class ApiIntragroup extends AtomsetIteratorAdapter implements
		AtomsetIteratorBasisDependent, AtomsetIteratorDirectable {

	/**
	 * Constructor makes iterator that must have basis specified and then be 
	 * reset() before iteration.
	 */
	public ApiIntragroup(Species[] species) {
		super(new ApiInnerVariable(
					new AtomIteratorBasis(),
					new AtomIteratorNbrCell(species[0],true)));
		pairIterator = (ApiInnerVariable)iterator;
		aiOuter = (AtomIteratorBasis)pairIterator.getOuterIterator();
		aiInner = (AtomIteratorSequencerList)pairIterator.getInnerIterator();
		aiInner.setNumToSkip(1);
	}

	/* (non-Javadoc)
	 * @see etomica.AtomsetIteratorBasisDependent#setDirective(etomica.IteratorDirective)
	 */
	public void setTarget(Atom[] targetAtoms) {
		aiOuter.setTarget(targetAtoms);
        oneTarget = targetAtoms[0] != null && (targetAtoms.length == 1 || targetAtoms[1] == null);
	}
	
	/**
	 * Puts iterator in a state to begin iteration.
	 */
	public void reset() {
	    doBoth = (direction == null) && oneTarget;
        upListNow = (direction == null) || (direction == IteratorDirective.UP); 
        
		if(upListNow || doBoth) {
			aiInner.setIterationDirection(IteratorDirective.UP);
			upListNow = true;
		}
		else {
			aiInner.setIterationDirection(IteratorDirective.DOWN);
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
	
	public void setDirection(IteratorDirective.Direction direction) {
        this.direction = direction;
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
		return atoms[0].seq.preceeds(atoms[1]) && list.contains(atoms[0]) && list.contains(atoms[1]);
	}

	/**
	 * Returns 1, indicating the the array given to the setBasis method
	 * should have only one element.
	 */
	public int basisSize() {
		return 1;
	}
	
    public AtomIteratorSequencerList getInnerIterator() {return aiInner;}
    
	protected final ApiInnerVariable pairIterator;//local, specifically typed copy
	private final AtomIteratorBasis aiOuter;//local, specifically typed copy
	private final AtomIteratorSequencerList aiInner;//local, specifically typed copy
	private boolean doBoth;//indicates if inner iteration should proceed UP and then DOWN from atom
    private boolean oneTarget;
    private IteratorDirective.Direction direction;
	private boolean upListNow;//indicates if inner iteration is currently in the UP direction
	private Atom basis;//atom most recently specified in setBasis; used by size() and contains(Atom[]) methods
}
