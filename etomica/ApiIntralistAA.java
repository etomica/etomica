package etomica;

/**
 * Loops through all molecule pairs given within a phase.
 *
 * @author David Kofke
 */
 
/* History 
 * 08/04/04 Created
 */

public class ApiIntralistAA implements AtomPairIterator, AtomPairListIterator {
    
    public ApiIntralistAA(Simulation sim) {
        pair = new AtomPair(sim.space);
        aiOuter = new AtomIteratorListSimple();
        aiInner = new AtomIteratorListSimple();
    }
    
	public void all(AtomSet basis, IteratorDirective id, final AtomsetActive action) {
		if(basis == null || !(action instanceof AtomPairActive)) return;
		switch(basis.nBody()) {
			case 1: all((Atom)basis, id, (AtomPairActive)action); break;
			case 2: all((AtomPair)basis, id, (AtomPairActive)action); break;
		}
	}
	public void all(Phase phase, IteratorDirective dummy, AtomPairActive action) {
		all(phase.speciesMaster().atomList, dummy, action);
	}
	
	public void all(AtomList basisList, IteratorDirective dummy, AtomPairActive action) {
		if(basisList == null || action == null) return;
		final AtomLinker header = basisList.header;
		AtomPair pair;
		for(AtomLinker e=header.next; e!=header; e=e.next) {
			if (e.atom != null) {
				for (AtomLinker f=e.next; f!=header; f=f.next) {
					if(f.atom != null) {
						pair = new AtomPair(e.atom, f.atom);
						action.actionPerformed(pair);
					}
				}
			}
		}
	}
	
	public void all(Atom basis, IteratorDirective dummy, AtomPairActive action) {
		throw new IllegalArgumentException("Error: ApiIntraphaseAA not defined for a pair basis");
/*		if(basis == null || action == null) return;
		AtomPairActive.OuterWrapper outerWrapper = action.outerWrapper();
		outerWrapper.setBoundary(basis.node.parentPhase().boundary());
		outerWrapper.aiInner = aiInner;
		outerWrapper.innerBasis = basis;
		outerWrapper.innerSkipFirst = true;
		aiOuter.all(basis, dummy, outerWrapper);*/
	}
	public void all(AtomPair basis, IteratorDirective id, AtomPairActive action) {
		throw new IllegalArgumentException("Error: ApiIntraphaseAA not defined for a pair basis");
	}

	public void setBasis(Phase phase) {
		setBasis(phase.speciesMaster().atomList);
	}

	public void setBasis(AtomList atomList) {
		aiOuter.setBasis(atomList);
		aiInner.setBasis(atomList);
		reset();
		pair.cPair.setBoundary(atomList.getFirst().node.parentPhase().boundary());
	}

	// this isn't an appropriate method for this class
	public void setBasis(Atom a1, Atom a2) {
        throw new IllegalArgumentException("Improper basis given to ApiIntraphaseAA");
   }
    
    /**
     * Returns the number of pairs capable of being given by this iterator.
     */
    public int size() {
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
//        pair.atom2 = aiInner.next();
//        pair.reset();
		pair.reset(pair.atom1,aiInner.next());//DAK 08/25/03  commented out two preceding lines
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
    
    protected boolean hasNext;
    protected boolean needUpdate1;
    
    protected final AtomIteratorListSimple aiOuter;
    protected final AtomIteratorListSimple aiInner;
    
    protected final IteratorDirective localDirective = new IteratorDirective(IteratorDirective.UP);
    protected final AtomPair pair;
    protected Atom atom1;
    
 //   private int count, outCount, inCount;//used in debugging to count the number of pairs given by the iterator
        
}  //end of class ApiIntraspeciesAA
    
