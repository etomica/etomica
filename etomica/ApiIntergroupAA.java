package etomica;

/**
 * Loops through all molecule pairs given between two groups.
 *
 * @author David Kofke
 */

/* History 
 * 08/25/03 (DAK) modified next method to invoke reset(Atom, Atom) on AtomPair.
 * Previously had pair.atom2 = iterator.next(); pair.reset().  Change made
 * because AtomPair's blank reset method does not call cPair.reset with atom
 * coordinate arguments.
 */


public final class ApiIntergroupAA implements AtomPairIterator {
    
    public ApiIntergroupAA(Simulation sim) {
        pair = new AtomPair(sim.space);
        aiOuter = sim.iteratorFactory.makeGroupIteratorSequential();
        aiInner = sim.iteratorFactory.makeIntergroupNbrIterator();
    }
    
	public void all(AtomSet basis, IteratorDirective id, final AtomsetActive action) {
		if(basis == null || !(action instanceof AtomPairActive)) return;
		switch(basis.nBody()) {
			case 1: all((Atom)basis, id, (AtomPairActive)action); break;
			case 2: all((AtomPair)basis, id, (AtomPairActive)action); break;
		}
	}
	public void all(Atom basis, IteratorDirective id, AtomPairActive action) {
		throw new IllegalArgumentException("Error: AtomPairIterator not defined for a single basis atom");
	}
	public void all(AtomPair basis, IteratorDirective dummy, AtomPairActive action) {
		if(basis == null || action == null) return;
		Atom group1 = basis.atom1();//assume group1 preceeds group2
		Atom group2 = basis.atom2();
		if(group1 == group2) throw new IllegalArgumentException("Improper basis given to ApiIntergroup1A: Basis atoms must be different");
		AtomPairActive.OuterWrapper outerWrapper = action.outerWrapper();
		outerWrapper.setBoundary(group1.node.parentPhase().boundary());
		outerWrapper.aiInner = aiInner;
		outerWrapper.innerBasis = group2;
		outerWrapper.innerSkipFirst = false;
		aiOuter.all(group1, dummy, outerWrapper);
	}

    public void setBasis(Atom a1, Atom a2) {
        if(a1 == a2)
            throw new IllegalArgumentException("Improper basis given to ApiInterGroupAA");
        group1 = a1;
        group2 = a2;
        aiOuter.setBasis(a1);
        aiInner.setBasis(a2);
		pair.cPair.setBoundary(group1.node.parentPhase().boundary());
    }
    
    /**
     * Returns the number of pairs capable of being given by this iterator.
     */
    public int size() {
        if(group1 == null || group2 == null) return 0;
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
        hasNext = false;
        if(group1 == null || group2 == null) return;

        aiOuter.reset();
        while(aiOuter.hasNext()) { //loop over iterator 1...
            pair.atom1 = aiOuter.next();
            aiInner.reset(localDirective.set(pair.atom1));
            if(aiInner.hasNext()) {
                hasNext = true;
                needUpdate1 = false;
                break;        //...until iterator 2 hasNext
            }
        }//end while
    }
        
    public AtomPair next() {
        //we use this update flag to indicate that atom1 in pair needs to be set to a new value.
        //it is not done directly in the while-loop because pair must first return with the old atom1 intact
        if(needUpdate1) {pair.atom1 = atom1; needUpdate1 = false;}  //aiOuter was advanced
//		pair.atom2 = aiInner.next();
//		pair.reset();
 	    pair.reset(pair.atom1,aiInner.next());//DAK 08/25/03  and commented out two preceding lines
        while(!aiInner.hasNext()) {     //Inner is done for this atom1, loop until it is prepared for next
            if(aiOuter.hasNext()) {     //Outer has another atom1...
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
    	throw new RuntimeException("ApiIntergroupAA.allPairs not impletmented");
//        act.outerWrapper.aiInner = aiInner;
//        aiOuter.reset();
//        aiOuter.allAtoms(outerWrapper);
//        hasNext = false;
    }

    private Atom group1, group2; 
    private boolean hasNext;
    private boolean needUpdate1;
    
    private final AtomIterator aiOuter;
    private final AtomIterator aiInner;
    
    //localdirective direction is set to BOTH to economize inner-loop iteration
    //with this setting, the inner iterator doesn't bother to check if its iterates
    //are up or down list of the reference atom (given by the outer loop).  Whatever
    //sets the basis for this pair iterator should ensure that the basis atoms
    //are compatible with the direction of iteration it wants
    private final IteratorDirective localDirective = new IteratorDirective(IteratorDirective.BOTH);
    private final AtomPair pair;
    private Atom atom1;
    
    
}  //end of class ApiIntergroupAA
    
