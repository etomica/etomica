package etomica;

/**
 * Basic class for iterating over pairs of atoms.
 * Pairs are iterated by collecting pairs yielded by two atom iterators.
 * Different types of pair iterators can be constructed with different choices
 * of the atom iterators.
 *
 * @author David Kofke
 */
public abstract class AtomPairIterator implements java.io.Serializable {
    
    private final AtomPair pair;
    /**
     * The iterators used to generate the sets of atoms
     */
    protected AtomIterator ai1, ai2;
    
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
    /**
     * Construct a pair iterator for use in the given phase.  Initial state is hasNext = false
     */
    public AtomPairIterator(Phase p) {
        pair = new AtomPair(p); 
        actionWrapper = new AtomPairAction.Wrapper(pair);
        hasNext = false;
    }
    /**
     * Construct a pair iterator for use in the given phase and which performs the given action.
     * Initial state is hasNext = false
     */
    public AtomPairIterator(Phase p, AtomPairAction.Wrapper wrap) {
        pair = new AtomPair(p); 
        actionWrapper = wrap;
        hasNext = false;
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

    /**
     * Resets the iterator, so that it is ready to go through all of its pairs.
     */
    public abstract void reset();

    /**
     * Resets the iterator so that it iterates over all pairs formed with the 
     * given atom.
     */
    public abstract void reset(Atom atom);

    /**
     * Resets iterator so that it iterates over all pairs formed from iterates
     * between and including the given atoms.
     */
    public abstract void reset(Atom atom1, Atom atom2);
        
        
    /**
        * Resets the first atom iterator with the given atom as an argument
        * Then resets the second iterator using the first non-null atom obtained from 
        * the first iterator.
        */
        //removed because of conflict with reset(IteratorDirective) and call to reset(null)
/*        public void reset(Atom a1) {
        if(a1==null) ai1.reset();
        else ai1.reset(a1);
        do {                  //advance up list of atoms until a pair is found
            if(ai1.hasNext()) {
                atom1 = ai1.next();
                ai2.reset(atom1);
            }  
            else {hasNext = false; return;}}   //end of list of atoms
        while(!ai2.hasNext());
        needUpdate1 = true;
        hasNext = true;
    }*/
            
    public final AtomPair next() {
        //why not do away with this line and set pair.atom1 in while loop?
///        if(needUpdate1) {pair.atom1 = atom1; needUpdate1 = false;}  //ai1 was advanced
        pair.atom2 = ai2.next();
        pair.reset();
        while(!ai2.hasNext()) {     //ai2 is done for this atom1, loop until it is prepared for next
            if(ai1.hasNext()) {     //ai1 has another atom1...
                atom1 = ai1.next();           //...get it
                if(ai2.reset(atom1) == atom1) //...reset ai2
                              ai2.next();     //...and advance if it's consequently set to return atom1
///                needUpdate1 = true;           //...flag update of pair.atom1 for next time
                pair.atom1 = atom1;
            }
            else {hasNext = false; break;} //ai1 has no more; all done with pairs
        }//end while
        return pair;
    }
    
    /**
     * This is called after ai1 and ai2 are defined and ai1 is reset, 
     * and it readies iterators to give first pair (or put hasNext = false if no pairs
     * are forthcoming).
     */
    protected final void setFirst() {
        while(ai1.hasNext()) { //loop over iterator 1...
            pair.atom1 = ai1.next();
            if(ai2.reset(pair.atom1) == pair.atom1) ai2.next(); //reset iterator 2 and advance if its first atom is atom1
            if(ai2.hasNext()) {
                hasNext = true;
                return;        //...until iterator 2 hasNext
            }
        }//end while
        hasNext = false;
    }//end setFirst
        
    /**
     * Performs the given action on all pairs returned by this iterator
     */
    public void allPairs(AtomPairAction act) {  
        reset();
        ai1.reset();  //this shouldn't be here, in general; need to consider it more carefully
        actionWrapper.pairAction = act;
        while(ai1.hasNext()) {
            pair.atom1 = ai1.next();
            ai2.reset(pair.atom1);
            ai2.allAtoms(actionWrapper);
        }
    }
}  //end of class AtomPairIterator
    
