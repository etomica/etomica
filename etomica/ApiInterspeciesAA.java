package etomica;

/**
 * Loops through all molecule pairs given within a single species.
 *
 * @author David Kofke
 */
public final class ApiInterspeciesAA implements AtomPairIterator {
    
    private Species species1, species2; 
    private boolean hasNext;
    private boolean needUpdate1;
    
    private final AtomIterator aiOuter = new AtomIteratorChildren();
    //this should be a neighbor iterator
    private final AtomIterator aiInner = new AtomIteratorChildren();
    
    private final IteratorDirective localDirective = new IteratorDirective();
    private final AtomPair pair;
    private Atom atom1;
    
    public ApiInterspeciesAA() {
        aiInner.setAsNeighbor(true);
        pair = new AtomPair(Simulation.instance.space);
    }
    
    public void setBasis(Atom a1, Atom a2) {
        if((a1 == a2) || !(a1 instanceof SpeciesAgent && a2 instanceof SpeciesAgent))
            throw new IllegalArgumentException("Improper basis given to ApiInterSpeciesAA");
        species1 = a1.node.parentSpecies();
        species2 = a2.node.parentSpecies();
        aiOuter.setBasis(a1);
        aiInner.setBasis(a2);
    }
    
    /**
     * Returns the number of pairs capable of being given by this iterator.
     */
    public int size() {
        if(species1 == null || species2 == null) return 0;
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
        if(species1 == null || species2 == null) {
            hasNext = false;
            return;
        }
        aiOuter.reset();
        setFirst();
    }
        
    /**
     * Resets the iterator so that it iterates over all pairs formed with the 
     * given atom in the most recently specified iterator directive (default UP is
     * if none previously specified.
     */
    public void reset(Atom atom) {  //need an exception
        System.out.println("Should not call reset(Atom) in ApiInterspeciesAA");
        System.exit(1);
    }
    
    public AtomPair next() {
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
     * This is called after aiOuter and aiInner are are reset, 
     * and it readies iterators to give first pair (or puts hasNext = false if no pairs
     * are forthcoming).
     */
    protected final void setFirst() {
        hasNext = false;
        while(aiOuter.hasNext()) { //loop over iterator 1...
            pair.atom1 = aiOuter.next();
            //keep this even though this is an interspecies iterator,
            //because neighbor iterator may reset based on atom1
            //right now, directive is ignored becuase aiInner is an instance
            //of AtomIteratorChildren
            aiInner.reset(localDirective.set(pair.atom1));
            if(aiInner.hasNext()) {
                hasNext = true;
                needUpdate1 = false;
                break;        //...until iterator 2 hasNext
            }
        }//end while
    }//end setFirst

    
    /**
     * Performs the given action on all pairs returned by this iterator.
     */
    public void allPairs(AtomPairAction act) {
        
    }
    
}  //end of class ApiInterspeciesAA
    
