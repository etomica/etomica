package etomica;

/**
 * Basic class for iterating over pairs of atoms.
 * Pairs are iterated by collecting atoms yielded by two atom iterators.
 * Different types of pair iterators can be constructed with different choices
 * of the atom iterators.
 *
 * @author David Kofke
 */
public class AtomPairIterator implements java.io.Serializable {
    
    private AtomPair pair; //want final, but Null inner class won't allow
    private IteratorDirective.Direction direction;
    private final IteratorDirective localDirective = new IteratorDirective();

    /**
     * The iterators used to generate the sets of atoms
     */
    protected /*final*/ AtomIterator ai1, ai2;
    protected AtomIterator aiInner, aiOuter;
    
    /**
     * A pair action wrapper used to enable the allPairs method
     */
    protected AtomPairAction.Wrapper actionWrapper;   // want final too //wrapper inner class defined below
    protected boolean hasNext;
    /**
     * Flag indicating whether atom1 of pair needs to be updated to point to the same atom that "atom1" in this class points to
     */
    private boolean needUpdate1; 
    private Atom atom1;
    
    //Used only for the NULL iterator, defined below.
    private AtomPairIterator() {
        hasNext = false;
        pair = null;
        ai1 = ai2 = null;
    }
    public AtomPairIterator(Phase p) {
        this(p, p.iteratorFactory().makeAtomIterator(), p.iteratorFactory().makeAtomIterator());
    }
    /**
     * Construct a pair iterator for the given phase, using the given atom iterators
     */
    public AtomPairIterator(Phase p, AtomIterator iter1, AtomIterator iter2) {
        pair = new AtomPair(p);
        actionWrapper = new AtomPairAction.Wrapper(pair);
        hasNext = false;
        ai1 = iter1;
        ai2 = iter2;
    }
    
    public final boolean hasNext() {return hasNext;}
        
    public void reset(IteratorDirective id) {
        direction = id.direction();
        switch(id.atomCount()) {
            case 0:  reset(); 
                     break;
            case 1:  reset(id.atom1()); 
                     break;
            case 2:  reset(id.atom1(), id.atom2()); 
                     break;
            default: hasNext = false; 
                     break;
        }
    }
    
    private final void setOuterInner(AtomIterator outIter, AtomIterator inIter) {
        aiOuter = outIter;
        aiInner = inIter;
        aiOuter.setAsNeighbor(false);
        aiInner.setAsNeighbor(true);
    }

    /**
     * Resets the iterator, so that it is ready to go through all of its pairs.
     */
    public void reset() {
        if(direction == IteratorDirective.BOTH) direction = IteratorDirective.UP;
        setOuterInner(ai1, ai2);
        aiOuter.reset(localDirective.set().set(direction));
        aiInner.reset(localDirective);
        setFirst();
    }
        
    /**
     * Resets the iterator so that it iterates over all pairs formed with the 
     * given atom.
     */
    public void reset(Atom atom) {
        if(ai1.contains(atom)) setOuterInner(ai1, ai2);
        else if(ai2.contains(atom)) setOuterInner(ai2, ai1);
        else {hasNext = false; return;}
        aiOuter.reset(localDirective.set(atom).set(IteratorDirective.NEITHER));
        aiInner.reset(localDirective.set().set(direction));//reset inner last to leave localDirective in proper state for inner-loop resets
        setFirst();
    }
    
    /**
     * Resets iterator so that it iterates over all pairs formed from iterates
     * between and including the given atoms.
     */
     public void reset(Atom atom1, Atom atom2) {
        ai1.setAsNeighbor(false);
        ai2.setAsNeighbor(false);
        if(ai1.reset(localDirective.set(atom1,atom2).set(direction)) != null) {
            setOuterInner(ai1, ai2);
        }
        else if(ai2.reset(localDirective) != null) {
            setOuterInner(ai2, ai1);
        }
        else {hasNext = false; return;}
        aiOuter.reset(localDirective);
        //maybe want to do something with setting inner considering atom2
        setFirst();
     }
        
    /**
     * This is called after aiOuter and aiInner are are reset, 
     * and it readies iterators to give first pair (or puts hasNext = false if no pairs
     * are forthcoming).
     */
    protected final void setFirst() {
 //       localDirective.set(IteratorDirective.UNCHANGED); //prepare for use in resetting inner loop here and in next
        hasNext = false;
        while(aiOuter.hasNext()) { //loop over iterator 1...
            pair.atom1 = aiOuter.next();
            aiInner.reset(localDirective.set(pair.atom1));
//            if(aiInner.reset(pair.atom1) == pair.atom1) aiInner.next(); //reset iterator 2 and advance if its first atom is atom1
            if(aiInner.hasNext()) {
                hasNext = true;
                needUpdate1 = false;
                break;        //...until iterator 2 hasNext
            }
        }//end while
    }//end setFirst
        
    public final AtomPair next() {
        //we use this update flag to indicate that atom1 in pair needs to be set to a new value.
        //it is not done directly in the while-loop because pair must first return with the old atom1 intact
        if(needUpdate1) {pair.atom1 = atom1; needUpdate1 = false;}  //aiOuter was advanced
        pair.atom2 = aiInner.next();
        pair.reset();
        while(!aiInner.hasNext()) {     //Inner is done for this atom1, loop until it is prepared for next
            if(aiOuter.hasNext()) {     //Outer has another atom1...
                atom1 = aiOuter.next();           //...get it
                aiInner.reset(localDirective.set(atom1)); //...reset Inner
                needUpdate1 = true;           //...flag update of pair.atom1 for next time
            }
            else {hasNext = false; break;} //Outer has no more; all done with pairs
        }//end while
        return pair;
    }

    /**
     * Performs the given action on all pairs returned by this iterator.
     */
     //not carefully checked
    public void allPairs(AtomPairAction act) {  
        reset();
//        ai1.reset();  //this shouldn't be here, in general; need to consider it more carefully
        actionWrapper.pairAction = act;
        while(ai1.hasNext()) {
            pair.atom1 = ai1.next();
            ai2.reset(localDirective.set(pair.atom1));
            ai2.allAtoms(actionWrapper);
        }
    }
    
    public static final AtomPairIterator NULL = new Null();
    private static final class Null extends AtomPairIterator {
        private Null() {super();}
        public void reset(IteratorDirective id) {}
    }
    
}  //end of class AtomPairIterator
    
