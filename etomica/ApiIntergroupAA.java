package etomica;

/**
 * Loops through all molecule pairs given between two groups.  Ignores
 * neighbor status and specification of direction.
 *
 * @author David Kofke
 */
 
public final class ApiIntergroupAA implements AtomPairIterator {
    
    public ApiIntergroupAA(Simulation sim) {
        pair = new AtomPair(sim.space);
        //inner and outer are arbitrary designations
        aiOuter = sim.iteratorFactory.makeGroupIteratorSimple();
        aiInner = sim.iteratorFactory.makeGroupIteratorSimple();
        outerWrapper = new AtomPairAction.OuterWrapper(pair, localDirective);
        outerWrapper.aiInner = aiInner;
    }
    
    public void setBasis(Atom a1, Atom a2) {
        aiOuter.setBasis(a1);
        aiInner.setBasis(a2);
    }
    
    /**
     * Returns the number of pairs capable of being given by this iterator.
     */
    public int size() {
        return aiOuter.size()*aiInner.size();
    }  
    
    public boolean hasNext() {return hasNext;}
    
    /**
     * Same as reset() -- iterator directive is ignored.
     */
    public void reset(IteratorDirective id) {
        reset();
    }
    
    /**
     * Resets the iterator, so that it is ready to go through all of its pairs.
     */
    public void reset() {
        aiOuter.reset();
        hasNext = aiOuter.hasNext();
        if(!hasNext) return;
        aiInner.reset();
        hasNext = aiInner.hasNext();
        if(!hasNext) return;
        pair.atom1 = aiOuter.next();
        needUpdate1 = false;
    }
        
    /**
     * Resets the iterator so that it iterates over all pairs formed with the 
     * given atom in the most recently specified iterator directive (default UP is
     * if none previously specified.
     */
    public void reset(Atom atom) {
        throw new IllegalArgumentException("Should not call reset(Atom) in ApiIntergroupAA");
    }
    
    public AtomPair next() {
        //we use this update flag to indicate that atom1 in pair needs to be set to a new value.
        //it is not done directly in the while-loop because pair must first return with the old atom1 intact
        if(needUpdate1) {pair.atom1 = atom1; needUpdate1 = false;}  //aiOuter was advanced
        pair.atom2 = aiInner.next();
        pair.reset();
        if(!aiInner.hasNext()) {
            if(aiOuter.hasNext()) {     //Outer has another atom1...
                atom1 = aiOuter.next();           //...get it
                aiInner.reset();  //don't pass directive since we don't care about neighbor iteration or direction
                needUpdate1 = true;           //...flag update of pair.atom1 for next time
            }
            else {hasNext = false;} //Outer has no more; all done with pairs
        }//end while
        return pair;
    }

    
    /**
     * Performs the given action on all pairs returned by this iterator.
     */
     //not carefully checked
    public void allPairs(AtomPairAction act) {
        //might be better to do loops explicitly, rather than 
        //via atom iterator
        outerWrapper.innerWrapper.pairAction = act;
        aiOuter.allAtoms(outerWrapper);
        hasNext = false;
    }
    
    private boolean hasNext;
    private boolean needUpdate1;
    
    //no neighbor iterator here (not species iterator)
    private final AtomIterator aiOuter;
    private final AtomIterator aiInner;
    
    private final AtomPairAction.OuterWrapper outerWrapper;
    private final IteratorDirective localDirective = new IteratorDirective();
    private final AtomPair pair;
    private Atom atom1;
    
}  //end of class ApiIntergroupAA
    
