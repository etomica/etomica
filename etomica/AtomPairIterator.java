package etomica;

/**
 * Basic class for iterating over pairs of atoms.
 * Pairs are iterated by collecting pairs yielded by two atom iterators.
 * Different types of pair iterators can be constructed with different choices
 * of the atom iterators.
 *
 * @author David Kofke
 */
public class AtomPairIterator implements java.io.Serializable {
    
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
    protected boolean needUpdate1; 
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
    public void reset() {
        reset(null);
    }
        
    //need to fill this in
    public void reset(IteratorDirective id) {
        
    }
        
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
    /**
        * Resets the first and second atom iterators using the first and second arguments, respectively.
        */
    public void reset(Atom a1, Atom a2) {
        ai1.reset(a1);
        ai2.reset(a2);
        pair.atom1 = ai1.next();
        needUpdate1 = false;
        hasNext = ai1.hasNext() && ai2.hasNext();
    }
            
    public AtomPair next() {
        if(needUpdate1) {pair.atom1 = atom1; needUpdate1 = false;}  //ai1 was advanced
        pair.atom2 = ai2.next();
        pair.reset();
        while(!ai2.hasNext()) {     //ai2 is done for this atom1, loop until it is prepared for next
            if(ai1.hasNext()) {     //ai1 has another atom1...
                atom1 = ai1.next();           //...get it
                if(ai2.reset(atom1) == atom1) //...reset ai2
                              ai2.next();     //...and advance if it's consequently set to return atom1
                needUpdate1 = true;           //...flag update of pair.atom1 for next time
            }
            else {hasNext = false; break;} //ai1 has no more; all done with pairs
        }//end while
        return pair;
    }
        
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
        
    // The following are convenience extensions of AtomPairIterator that
    // handle some common iteration needs
        
    /**
        * Iterator for all atom pairs in a phase
        * Default is to do inter and intra pairs; this may be overridden using reset method to do
        * only intermolecular pairs
        * Uses atom iterator and atomPair iterator given by the phase.iteratorFactory class.
        */
    public static final class All extends AtomPairIterator {
        public All(Phase p) {
            super(p);
            ai1 = p.iteratorFactory().makeAtomIteratorUp();
            ai2 = p.iteratorFactory().makeAtomIteratorUpNeighbor();
            this.reset();
        }
    }
         
    /**
    * Iterates over pairs formed by given atom and all atoms from other molecules above it in list
    * If given atom is not in phase, it is considered the last atom, and no iterations are performed
    */
    public static final class Up extends AtomPairIterator {
        public Up(Phase p) {
            super(p);
            ai1 = new AtomIterator.Singlet();
            ai2 = p.iteratorFactory().makeAtomIteratorUpNeighbor();
            this.reset();
        }
        public Up(Phase p, Atom a) {
            super(p);
            ai1 = new AtomIterator.Singlet(a);
            ai2 = p.iteratorFactory().makeAtomIteratorUpNeighbor();
            this.reset();
        }
    }
         
    /**
    * Iterates over pairs formed by given atom and all atoms from other molecules below it in list
    * If given atom is not in phase, it is considered the last atom, and iterations are performed
    * over pairs formed from it and all atoms in phase
    */
    public static final class Down extends AtomPairIterator {
        public Down(Phase p) {
            super(p);
            ai1 = new AtomIterator.Singlet();
            ai2 = p.iteratorFactory().makeAtomIteratorDownNeighbor();
            this.reset();
        }
        public Down(Phase p, Atom a) {
            super(p);
            ai1 = new AtomIterator.Singlet(a);
            ai2 = p.iteratorFactory().makeAtomIteratorDownNeighbor();
            this.reset();
        }
    }        
}  //end of class AtomPairIterator
    
