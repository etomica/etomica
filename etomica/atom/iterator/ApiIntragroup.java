/*
 * History
 * Created on Aug 30, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.action.AtomsetCount;
import etomica.action.AtomsetDetect;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;

/**
 * Returns iterates from from the childList of a single basis atom.  Behavior is set
 * via iterator directive: if no atom is specified there, all pairs formed from the childList
 * are given; otherwise, if an atom is specified there, pairs will be formed from the
 * childList atoms with the basis' child from which the target atom is descended.  
 */
public final class ApiIntragroup implements AtomPairIterator,
		AtomsetIteratorBasisDependent, AtomsetIteratorDirectable {

	/**
	 * Constructor makes iterator that must have basis specified and then be 
	 * reset() before iteration.
	 */
	public ApiIntragroup() {
		this(new AtomIteratorSequence(IteratorDirective.UP, 1),
                new AtomIteratorSequence(IteratorDirective.DOWN, 1));
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
    
    public void setTarget(AtomSet targetAtoms) {
        if(targetAtoms == null) throw new NullPointerException("Cannot set target to null; use AtomSet.NULL");
        aiOuter.setTarget(targetAtoms);
        oneTarget = targetAtoms.count() != 0 && targetAtoms.getAtom(0) != null 
                        && (targetAtoms.count() == 1 || targetAtoms.getAtom(1) == null);
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
        //upList if no target given (then do all pairs) or if specified by direction
        upListNow = (!oneTarget || direction != IteratorDirective.DOWN);
        //dnList only if one target and not explicitly directed up
        doGoDown = (oneTarget && direction != IteratorDirective.UP);
        
        if (upListNow) {
            apiUp.reset();
        }
        else {//if upListNow is false, doGoDown must be true
            apiUp.unset();
            apiDown.reset();
        }
	}
    
    public void unset() {
        apiDown.unset();
        apiUp.unset();
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
    
    public boolean hasNext() {
        return upListNow ? apiUp.hasNext() : apiDown.hasNext();
    }

    public AtomSet next() {
        return nextPair();
    }
    
    public AtomSet peek() {
        return upListNow ? apiUp.peek() : apiDown.peek();
    }
    
    public AtomPair nextPair() {
        if (upListNow) {
            AtomPair next = apiUp.nextPair();
            if (!apiUp.hasNext() && doGoDown) {
                upListNow = false;
                apiDown.reset();
            }
            return next;
        }
        return apiDown.nextPair();
    }
    
    /**
     * Performs action on all iterates for iterator as currently
     * conditioned (basis, target, direction).  Clobbers
     * iteration state.
     */
    public void allAtoms(AtomsetAction action) {
        if (!oneTarget || direction != IteratorDirective.DOWN) {
            apiUp.allAtoms(action);
        }
        if (oneTarget && direction != IteratorDirective.UP) {
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
    
	private final AtomsetIteratorBasisDependent aiOuter;//local, specifically typed copy
	private boolean oneTarget;
    private IteratorDirective.Direction direction;
    private AtomsetCount counter;
    private AtomsetDetect detector;
    private final ApiInnerVariable apiUp, apiDown;
    private boolean doGoDown, upListNow;
}
