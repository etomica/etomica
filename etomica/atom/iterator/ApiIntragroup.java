/*
 * History
 * Created on Aug 30, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.IteratorDirective;
import etomica.action.AtomsetAction;
import etomica.action.AtomsetCount;
import etomica.action.AtomsetDetect;
import etomica.atom.AtomIterator;
import etomica.atom.AtomSet;

/**
 * Returns iterates from from the childList of a single basis atom.  Behavior is set
 * via iterator directive: if no atom is specified there, all pairs formed from the childList
 * are given; otherwise, if an atom is specified there, pairs will be formed from the
 * childList atoms with the basis' child from which the target atom is descended.  
 */
public final class ApiIntragroup extends AtomPairIteratorAdapter implements
		AtomsetIteratorBasisDependent, AtomsetIteratorDirectable, ApiComposite {

	/**
	 * Constructor makes iterator that must have basis specified and then be 
	 * reset() before iteration.
	 */
	public ApiIntragroup() {
		this(new ApiInnerVariable(
					new AtomIteratorBasis(),
					new AtomIteratorSequenceDirectable()));
        ((AtomIteratorSequenceDirectable)aiInner).setNumToSkip(1);
    }
    
    public ApiIntragroup(ApiComposite pairIterator) {
        super(pairIterator);
		aiOuter = (AtomIteratorBasis)pairIterator.getOuterIterator();
		aiInner = (AtomsetIteratorDirectable)pairIterator.getInnerIterator();
	}

	public void setTarget(AtomSet targetAtoms) {
        if(targetAtoms == null) throw new NullPointerException("Cannot set target to null; use AtomSet.NULL");
		aiOuter.setTarget(targetAtoms);
        oneTarget = targetAtoms.count() != 0 && targetAtoms.getAtom(0) != null && (targetAtoms.count() == 1 || targetAtoms.getAtom(1) == null);
	}
    
    public boolean haveTarget(AtomSet targetAtoms) {
        for(int i=targetAtoms.count()-1; i>=0; i--) {
            if(!aiOuter.haveTarget(targetAtoms.getAtom(i))) return false;
        }
        return true;
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
	public void setBasis(AtomSet atoms) {
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
        if(counter == null) {
            counter = new AtomsetCount();
        } else {
            counter.reset();
        }
        allAtoms(counter);
        return counter.callCount();
	}
    
    public boolean contains(AtomSet atoms) {
        if(atoms == null || atoms.count() != 2) return false;
        if(detector == null) {
            detector = new AtomsetDetect(atoms);
        } else {
            detector.setAtoms(atoms);
            detector.reset();
        }
        allAtoms(detector);
        return detector.detectedAtom();
    }
    
    /**
     * Performs action on all iterates for iterator as currently
     * conditioned (basis, target, direction).  Clobbers
     * iteration state.
     */
    public void allAtoms(AtomsetAction action) {
        //iterate in the prescribed direction from the target
        if(oneTarget) aiInner.setDirection(direction);
        
        //no target given -- iterate over all pairs
        else aiInner.setDirection(IteratorDirective.UP);
        
        super.allAtoms(action);
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
	private final AtomsetIteratorDirectable aiInner;//local, specifically typed copy
	private boolean oneTarget;
    private IteratorDirective.Direction direction;
    private AtomsetCount counter;
    private AtomsetDetect detector;
}
