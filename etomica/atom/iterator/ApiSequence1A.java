/*
 * History
 * Created on Aug 30, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.action.AtomsetAction;
import etomica.atom.Atom;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;

/**
 * Returns iterates from a list.  
 */
public class ApiSequence1A implements AtomPairIterator, AtomsetIteratorDirectable {

    /**
     * Constructor makes iterator that must have list specified and then be 
     * reset() before iteration.
     */
    public ApiSequence1A() {
        this(new AtomIteratorSequence(IteratorDirective.UP, 1),
                new AtomIteratorSequence(IteratorDirective.DOWN, 1));
    }
    
    public ApiSequence1A(AtomIteratorAtomDependent aiInnerUp, AtomIteratorAtomDependent aiInnerDn) {
        this.aiOuter = new AtomIteratorSinglet();
        apiUp = new ApiInnerVariable(aiOuter, aiInnerUp, false);
        apiDown = new ApiInnerVariable(aiOuter, aiInnerDn, true);
    }

    public int nBody() {
        return 2;
    }
    
    public void setAtom(Atom atom) {
        aiOuter.setAtom(atom);
    }
    
    /**
     * Puts iterator in a state to begin iteration.
     */
    public void reset() {
        //upList if specified by direction
        upListNow = (direction != IteratorDirective.DOWN);
        //dnList only if not explicitly directed up
        doGoDown = (direction != IteratorDirective.UP);
        
        if (upListNow) {
            apiUp.reset();
            if (!apiUp.hasNext()) {
                upListNow = false;
            }
        }
        if (!upListNow && doGoDown) {
            apiUp.unset();
            apiDown.reset();
        }
    }
    
    public void unset() {
        apiDown.unset();
        apiUp.unset();
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
        if (aiOuter.getAtom() == null) {
            return 0;
        }
        int count = 0;
        if (direction != IteratorDirective.DOWN) {
            count += apiUp.size();
        }
        if (direction != IteratorDirective.UP) {
            count += apiDown.size();
        }
        return count;
    }
    
    public boolean contains(AtomSet atoms) {
        if(atoms == null || atoms.count() != 2) return false;
        if (direction != IteratorDirective.DOWN
                && apiUp.contains(atoms)) {
            return true;
        }
        if (direction != IteratorDirective.UP
                && apiDown.contains(atoms)) {
            return true;
        }
        return false;
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
        if (direction != IteratorDirective.DOWN) {
            apiUp.allAtoms(action);
        }
        if (direction != IteratorDirective.UP) {
            apiDown.allAtoms(action);
        }
    }
	
    private final AtomIteratorSinglet aiOuter;
    private IteratorDirective.Direction direction;
    private final ApiInnerVariable apiUp, apiDown;
    private boolean doGoDown, upListNow;
}
