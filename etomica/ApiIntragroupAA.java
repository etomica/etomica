package etomica;

/**
 * Loops through all molecule pairs given within a single group.
 *
 * @author David Kofke
 */
 
 //should put AA iterators into a single class, configured at construction
 //with appropriate inner integrators
 //this class and ApiInterspeciesAA are nearly identical
public final class ApiIntragroupAA extends AtomPairIterator {
    
    public ApiIntragroupAA(Simulation sim) {
        pair = new AtomPair(sim.space);
        aiOuter = sim.iteratorFactory.makeGroupIteratorSequential();
        aiInner = sim.iteratorFactory.makeIntragroupNbrIterator();
    }
    
	public void all(Atom basis, IteratorDirective dummy, AtomPairActive action) {
		if(basis == null || action == null) return;
		AtomPairActive.OuterWrapper outerWrapper = action.outerWrapper();
		outerWrapper.aiInner = aiInner;
		outerWrapper.innerBasis = basis;
		outerWrapper.innerSkipFirst = true;
		aiOuter.all(basis, dummy, outerWrapper);
	}
	public void all(AtomPair basis, IteratorDirective id, AtomPairActive action) {
		throw new IllegalArgumentException("Error: AtomPairIterator not defined for a pair basis");
	}

   public void setBasis(Atom a1, Atom a2) {
        if(a1 != a2)
            throw new IllegalArgumentException("Improper basis given to ApiIntraSpeciesAA");
        group = a1;
        aiOuter.setBasis(a1);
        aiInner.setBasis(a1);
    }
    
    /**
     * Returns the number of pairs capable of being given by this iterator.
     */
    public int size() {
        if(group == null) return 0;
        int n = aiOuter.size();
        return n*(n-1)/2;
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
        hasNext = false;
        if(group == null) return;

        aiOuter.reset();
        while(aiOuter.hasNext()) { //loop over iterator 1...
            pair.atom1 = aiOuter.next();
//            outCount++;
            aiInner.reset(localDirective.set(pair.atom1));
            if(aiInner.hasNext()) {
                hasNext = true;
                needUpdate1 = false;
                break;        //...until iterator 2 hasNext
            }
        }//end while
        
 //        System.out.println(count + "  " + outCount);
 //       count = outCount = 0;
    }
        
    public AtomPair next() {
  //      count++;
        //we use this update flag to indicate that atom1 in pair needs to be set to a new value.
        //it is not done directly in the while-loop because pair must first return with the old atom1 intact
        if(needUpdate1) {pair.atom1 = atom1; needUpdate1 = false;}  //aiOuter was advanced
        pair.atom2 = aiInner.next();
        pair.reset();
        while(!aiInner.hasNext()) {     //Inner is done for this atom1, loop until it is prepared for next
            if(aiOuter.hasNext()) {     //Outer has another atom1...
  //          outCount++;
                atom1 = aiOuter.next();           //...get it
                aiInner.reset(localDirective.set(atom1)); //...reset Inner (don't advance because it is inter-group)
                needUpdate1 = true;           //...flag update of pair.atom1 for next time
            }
            else {hasNext = false; break;} //Outer has no more; all done with pairs
        }//end while
        return pair;
    }

    /**
     * Performs the given action on all pairs returned by this iterator.
     */
    public void allPairs(AtomPairAction act) { 
    	throw new etomica.exception.MethodNotImplementedException();
//        act.outerWrapper.aiInner = aiInner;
//        aiOuter.reset();
//        aiOuter.allAtoms(act.outerWrapper);
//        hasNext = false;
    }
    
    private Atom group; 
    private boolean hasNext;
    private boolean needUpdate1;
    
    private final AtomIterator aiOuter;
    private final AtomIterator aiInner;
    
    private final IteratorDirective localDirective = new IteratorDirective(IteratorDirective.UP);
    private final AtomPair pair;
    private Atom atom1;
    
 //   private int count, outCount, inCount;//used in debugging to count the number of pairs given by the iterator
        
}  //end of class ApiIntraspeciesAA
    
