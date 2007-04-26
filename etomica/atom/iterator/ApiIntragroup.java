/*
 * History
 * Created on Aug 30, 2004 by kofke
 */
package etomica.atom.iterator;

import java.io.Serializable;

import etomica.action.AtomsetAction;
import etomica.action.AtomsetCount;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;

/**
 * Returns iterates from the childList of a single basis atom.  Behavior is set
 * via iterator directive: if no atom is specified there, all pairs formed from the childList
 * are given; otherwise, if an atom is specified there, pairs will be formed from the
 * childList atoms with the basis' child from which the target atom is descended.  
 */
public final class ApiIntragroup implements AtomPairIterator,
		AtomsetIteratorBasisDependent, AtomsetIteratorDirectable,
        Serializable {

    /**
	 * Constructor makes iterator that must have basis specified and then be 
	 * reset() before iteration.
	 */
	public ApiIntragroup() {
		this(new AtomIteratorArrayList(IteratorDirective.Direction.UP, 1),
                new AtomIteratorArrayList(IteratorDirective.Direction.DOWN, 1));
    }
    
    public ApiIntragroup(AtomIteratorAtomDependent aiInnerUp, AtomIteratorAtomDependent aiInnerDn) {
        this(new AtomIteratorBasis(), aiInnerUp, aiInnerDn);
    }
    
    public ApiIntragroup(AtomsetIteratorBasisDependent aiOuter, 
            AtomIteratorAtomDependent aiInnerUp, AtomIteratorAtomDependent aiInnerDn) {
        this.aiOuter = aiOuter;
        apiUp = new ApiInnerVariable((AtomIterator)aiOuter, aiInnerUp, false);
        apiDown = new ApiInnerVariable((AtomIterator)aiOuter, aiInnerDn, true);
	}

    public int nBody() {
        return 2;
    }
    
    public void setTarget(IAtom newTargetAtom) {
        targetAtom = newTargetAtom;
        aiOuter.setTarget(targetAtom);
	}
    
    public boolean haveTarget(IAtom target) {
        return aiOuter.haveTarget(target);
    }
	
	/**
	 * Puts iterator in a state to begin iteration.
	 */
	public void reset() {
        //upList if no target given (then do all pairs) or if specified by direction
        upListNow = (targetAtom == null || direction != IteratorDirective.Direction.DOWN);
        //dnList only if one target and not explicitly directed up
        doGoDown = (targetAtom != null && direction != IteratorDirective.Direction.UP);
        
        if (upListNow) {
            apiUp.reset();
        }
        else {
            apiDown.reset();
        }
	}
    
    public void unset() {
        upListNow = false;
        apiDown.unset();
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
    
    public AtomSet next() {
        return nextPair();
    }
    
    public AtomPair nextPair() {
        if (upListNow) {
            AtomPair next = apiUp.nextPair();
            if (next != null && (next.atom0 == null || next.atom1 == null)) {
                throw new RuntimeException("oops "+next);
            }
            if (next != null || !doGoDown) {
                return next;
            }
            upListNow = false;
            apiDown.reset();
        }
        return apiDown.nextPair();
    }
    
    /**
     * Performs action on all iterates for iterator as currently
     * conditioned (basis, target, direction).  Clobbers
     * iteration state.
     */
    public void allAtoms(AtomsetAction action) {
        if (targetAtom == null || direction != IteratorDirective.Direction.DOWN) {
            apiUp.allAtoms(action);
        }
        if (targetAtom != null && direction != IteratorDirective.Direction.UP) {
            apiDown.allAtoms(action);
        }
    }
	
	/**
	 * Returns 1, indicating the the array given to the setBasis method
	 * should have only one element.
	 */
	public int basisSize() {
		return 1;
	}
    
    private static final long serialVersionUID = 2L;
	private final AtomsetIteratorBasisDependent aiOuter;//local, specifically typed copy
	private IAtom targetAtom;
    private IteratorDirective.Direction direction;
    private AtomsetCount counter;
    private final ApiInnerVariable apiUp, apiDown;
    private boolean doGoDown, upListNow;
}
