package etomica;

/**
 * Forms pairs from 1 atom of one group with all the atoms from another group.
 *
 * @author David Kofke
 */
public class ApiIntergroup1A implements AtomPairIterator {
    
    private final AtomIteratorChildren aiInner = new AtomIteratorChildren();
    private final AtomIteratorChildren aiOuterA = new AtomIteratorChildren();
    private final AtomIteratorSinglet  aiOuter1 = new AtomIteratorSinglet();
    private AtomIterator aiOuter;
    private boolean hasNext;
    
    public void setBasis(Atom a1, Atom a2) {
        aiOuterA.setBasis(a1);
        aiInner.setBasis(a2);
    }
    
    /**
     * Returns the number of pairs capable of being given by this iterator
     * (that is, if no restrictions are specified in an iteratorDirective).
     */
    public int size() {return aiOuter.size()*aiInner.size();}       
    
    public boolean hasNext() {return hasNext;}
    
    public void reset(IteratorDirective id) {
        if(id.atom1() == null) {
            aiOuter = aiOuterA;
        } else {
            aiOuter = aiOuter1;
        }
    }
    
    /**
     * Resets the iterator, so that it is ready to go through all of its pairs.
     */
    public void reset() {
        aiOuter.reset();
        aiInner.reset();
    }
        
    /**
     * Resets the iterator so that it iterates over all pairs formed with the 
     * given atom.
     */
    public void reset(Atom atom);
    
        
    public AtomPair next();

    /**
     * Performs the given action on all pairs returned by this iterator.
     */
    public void allPairs(AtomPairAction act);
    
}  //end of class AtomPairIterator
    
