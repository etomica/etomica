package etomica;

/**
 * Pair iterator synthesized from two atom iterators.  Pairs are formed from
 * the atoms yielded by the two atom iterators.
 * Different types of pair iterators can be constructed with different choices
 * of the atom iterators.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 8/4/02 (DAK) special-purpose modification to setBasis method to in attempt to work with PistonCylinder
  * 08/29/03 (DAK) using reset2 method
  */
  
public final class ApiGeneral implements AtomPairIterator {
    
    private final AtomPair pair;
    private IteratorDirective.Direction direction;
    private final IteratorDirective localDirective = new IteratorDirective();

    /**
     * The iterators used to generate the sets of atoms
     */
    private final AtomIterator aiInner, aiOuter;
    
    protected boolean hasNext;
    /**
     * Flag indicating whether atom1 of pair needs to be updated to point to the same atom that "atom1" in this class points to
     */
    private boolean needUpdate1; 
    private Atom atom1;
    
    /**
     * Construct a pair iterator using the given atom iterators
     */
    public ApiGeneral(Space s, AtomIterator aiOuter, AtomIterator aiInner) {
        pair = new AtomPair(s);
        hasNext = false;
        this.aiOuter = aiOuter;
        this.aiInner = aiInner;
        direction = IteratorDirective.UP;
    }
 
	public void all(AtomSet basis, IteratorDirective id, final AtomSetActive action) {
		if(basis == null || !(action instanceof AtomPairActive)) return;
		switch(basis.nBody()) {
			case 1: all((Atom)basis, id, (AtomPairActive)action); break;
			case 2: all((AtomPair)basis, id, (AtomPairActive)action); break;
		}
	}
	public void all(Atom basis, IteratorDirective dummy, AtomPairActive action) {
		throw new RuntimeException("Method all not implemented in AtomIteratorCompound");
		//in progress
	}
	public void all(AtomPair basis, IteratorDirective dummy, AtomPairActive action) {
		if(basis == null || action == null) return;
		Atom group1 = basis.atom1();//assume group1 preceeds group2
		Atom group2 = basis.atom2();
		//ugly workaround for PistonCylinder
		if(group1.node.isLeaf() && aiInner instanceof AtomIteratorSinglet) {
			group1 = basis.atom2();
			group2 = basis.atom1();
		}
		AtomPairActive.OuterWrapper outerWrapper = action.outerWrapper();
		outerWrapper.setBoundary(group1.node.parentPhase().boundary());
		outerWrapper.aiInner = aiInner;
		outerWrapper.innerBasis = group2;
		aiOuter.all(group1, dummy, outerWrapper);
	}

   
    public void setBasis(Atom a1, Atom a2) {
		//ugly workaround for PistonCylinder
		if(a1.node.isLeaf() && aiInner instanceof AtomIteratorSinglet) {
			aiOuter.setBasis(a2);
			aiInner.setBasis(a1);
		} else {
			//originally was just this, without if/else block
			aiOuter.setBasis(a1);
			aiInner.setBasis(a2);
		}
		pair.cPair.setBoundary(a1.node.parentPhase().boundary());
    }
    
    public AtomIterator aiOuter() {return aiOuter;}
    public AtomIterator aiInner() {return aiInner;}
    
    /**
     * Returns the number of pairs capable of being given by this iterator
     * (that is, if no restrictions are specified in an iteratorDirective).
     */
    public int size() {
        throw new RuntimeException("ApiGeneral.size should be checked before using");
      /*  int n1 = aiOuter.size();
        int n2 = aiInner.size();
        return (aiOuter.getBasis() == aiInner.getBasis()) ?
            n1*(n2-1)/2 : n1*n2;*/
    }        
    
    public final boolean hasNext() {return hasNext;}
        
    public void reset(IteratorDirective id) {
        direction = id.direction();
        switch(id.atomCount()) {
            case 0:  reset(); 
                     break;
            case 1:  reset(id.atom1()); 
                     break;
            default: hasNext = false; 
                     break;
        }
    }

    /**
     * Resets the iterator, so that it is ready to go through all of its pairs.
     */
    public void reset() {
        aiOuter.reset();
        setFirst();
    }
        
    /**
     * Resets the iterator so that it iterates over all pairs formed with the 
     * given atom in the most recently specified iterator directive (default UP is
     * if none previously specified).
     */
    public void reset(Atom atom) {
        aiOuter.reset(localDirective.set(atom).set(IteratorDirective.NEITHER));
        localDirective.set(direction);
//        aiInner.reset(localDirective.set().set(direction));//reset inner last to leave localDirective in proper state for inner-loop resets
        setFirst();
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
            Atom nextInner = aiInner.reset(localDirective.set(pair.atom1));
            if(nextInner == pair.atom1) aiInner.next();//move to atom after atom1
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
//        pair.atom2 = aiInner.next();
        pair.reset2(aiInner.next());//08/29/03 deleted above line and added this
        while(!aiInner.hasNext()) {     //Inner is done for this atom1, loop until it is prepared for next
            if(aiOuter.hasNext()) {     //Outer has another atom1...
                atom1 = aiOuter.next();           //...get it
                Atom nextInner = aiInner.reset(localDirective.set(atom1)); //...reset Inner
                if(nextInner == atom1) aiInner.next();
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
        throw new RuntimeException("ApiGeneral.allPairs should be checked before using");
        /*outerWrapper.innerWrapper.pairAction = act;
        outerWrapper.aiInner = aiInner;
        aiOuter.allAtoms(outerWrapper);
        hasNext = false;
        */
    /*    while(ai1.hasNext()) {
            pair.atom1 = ai1.next();
            ai2.reset(localDirective.set(pair.atom1));
            ai2.allAtoms(actionWrapper);
        }*/
    }
    
/*    public static void main(String[] args) throws java.io.IOException {
        
        java.io.BufferedReader in = new java.io.BufferedReader(new java.io.InputStreamReader(System.in));
        Simulation.instance = new Simulation();
        Phase phase = new Phase();
        Phase phase2 = new Phase();
        Species species = new SpeciesSpheres();
        species.setNMolecules(8);
        Simulation.instance.elementCoordinator.go();

        //Phase phase = TestAll.setupTestPhase(4);
        Atom  atom = ((AtomGroup)phase.firstSpecies().getAtom(2)).firstChild();
     //   Atom  atom = phase.firstSpecies().getAtom(2);
        AtomPairIterator iterator = new AtomPairIterator(phase);
        IteratorDirective id = new IteratorDirective();
        String line;
        
        System.out.println("reset()"); iterator.reset();
        while(iterator.hasNext()) {AtomPair pair = iterator.next();System.out.println(pair.atom1().toString()+ "  " + pair.atom2().toString());}  line = in.readLine();

        System.out.println("reset(DOWN)"); iterator.reset(id.set(IteratorDirective.DOWN));
        while(iterator.hasNext()) {AtomPair pair = iterator.next();System.out.println(pair.atom1().toString()+ "  " + pair.atom2().toString());}  line = in.readLine();

        System.out.println("reset(NEITHER)"); iterator.reset(id.set(IteratorDirective.NEITHER));
        while(iterator.hasNext()) {AtomPair pair = iterator.next();System.out.println(pair.atom1().toString()+ "  " + pair.atom2().toString());}  line = in.readLine();
 
        System.out.println("reset(BOTH)"); iterator.reset(id.set(IteratorDirective.BOTH));
        while(iterator.hasNext()) {AtomPair pair = iterator.next();System.out.println(pair.atom1().toString()+ "  " + pair.atom2().toString());}  line = in.readLine();
 
        System.out.println("reset(atom,UP)"); iterator.reset(id.set(atom).set(IteratorDirective.UP));
        while(iterator.hasNext()) {AtomPair pair = iterator.next();System.out.println(pair.atom1().toString()+ "  " + pair.atom2().toString());}  line = in.readLine();
 
        System.out.println("reset(atom,DOWN)"); iterator.reset(id.set(atom).set(IteratorDirective.DOWN));
        while(iterator.hasNext()) {AtomPair pair = iterator.next();System.out.println(pair.atom1().toString()+ "  " + pair.atom2().toString());}  line = in.readLine();
 
        System.out.println("reset(atom,NEITHER)"); iterator.reset(id.set(atom).set(IteratorDirective.NEITHER));
        while(iterator.hasNext()) {AtomPair pair = iterator.next();System.out.println(pair.atom1().toString()+ "  " + pair.atom2().toString());}  line = in.readLine();
 
        System.out.println("reset(atom,BOTH)"); iterator.reset(id.set(atom).set(IteratorDirective.BOTH));
        while(iterator.hasNext()) {AtomPair pair = iterator.next();System.out.println(pair.atom1().toString()+ "  " + pair.atom2().toString());}  line = in.readLine();

        atom = phase.lastAtom();
        
        System.out.println("reset(lastatom,UP)"); iterator.reset(id.set(atom).set(IteratorDirective.UP));
        while(iterator.hasNext()) {AtomPair pair = iterator.next();System.out.println(pair.atom1().toString()+ "  " + pair.atom2().toString());}  line = in.readLine();
 
        System.out.println("reset(lastatom,DOWN)"); iterator.reset(id.set(atom).set(IteratorDirective.DOWN));
        while(iterator.hasNext()) {AtomPair pair = iterator.next();System.out.println(pair.atom1().toString()+ "  " + pair.atom2().toString());}  line = in.readLine();
 
        System.out.println("reset(lastatom,NEITHER)"); iterator.reset(id.set(atom).set(IteratorDirective.NEITHER));
        while(iterator.hasNext()) {AtomPair pair = iterator.next();System.out.println(pair.atom1().toString()+ "  " + pair.atom2().toString());}  line = in.readLine();
 
        System.out.println("reset(lastatom,BOTH)"); iterator.reset(id.set(atom).set(IteratorDirective.BOTH));
        while(iterator.hasNext()) {AtomPair pair = iterator.next();System.out.println(pair.atom1().toString()+ "  " + pair.atom2().toString());}  line = in.readLine();

        Species species2 = new SpeciesSpheres();
        species.setNMolecules(4);
        species2.setNMolecules(5);
        Simulation.instance.elementCoordinator.go();
        
        System.out.println("second species added");
        iterator = new AtomPairIterator(phase,
            species.getAgent(phase).makeLeafAtomIterator(), 
            species2.getAgent(phase).makeLeafAtomIterator());
 
        atom = ((AtomGroup)species2.getAgent(phase).getAtom(3)).firstChild();
        
        System.out.println("reset(atom,UP)"); iterator.reset(id.set(atom).set(IteratorDirective.UP));
        while(iterator.hasNext()) {AtomPair pair = iterator.next();System.out.println(pair.atom1().toString()+ "  " + pair.atom2().toString());}  line = in.readLine();
 
        System.out.println("reset(atom,DOWN)"); iterator.reset(id.set(atom).set(IteratorDirective.DOWN));
        while(iterator.hasNext()) {AtomPair pair = iterator.next();System.out.println(pair.atom1().toString()+ "  " + pair.atom2().toString());}  line = in.readLine();
 
        System.out.println("reset(atom,NEITHER)"); iterator.reset(id.set(atom).set(IteratorDirective.NEITHER));
        while(iterator.hasNext()) {AtomPair pair = iterator.next();System.out.println(pair.atom1().toString()+ "  " + pair.atom2().toString());}  line = in.readLine();
 
        System.out.println("reset(atom,BOTH)"); iterator.reset(id.set(atom).set(IteratorDirective.BOTH));
        while(iterator.hasNext()) {AtomPair pair = iterator.next();System.out.println(pair.atom1().toString()+ "  " + pair.atom2().toString());}  line = in.readLine();

        atom = ((AtomGroup)species.getAgent(phase).getAtom(3)).firstChild();
        
        System.out.println("reset(atom,UP)"); iterator.reset(id.set(atom).set(IteratorDirective.UP));
        while(iterator.hasNext()) {AtomPair pair = iterator.next();System.out.println(pair.atom1().toString()+ "  " + pair.atom2().toString());}  line = in.readLine();
 
        System.out.println("reset(atom,DOWN)"); iterator.reset(id.set(atom).set(IteratorDirective.DOWN));
        while(iterator.hasNext()) {AtomPair pair = iterator.next();System.out.println(pair.atom1().toString()+ "  " + pair.atom2().toString());}  line = in.readLine();
 
        System.out.println("reset(atom,NEITHER)"); iterator.reset(id.set(atom).set(IteratorDirective.NEITHER));
        while(iterator.hasNext()) {AtomPair pair = iterator.next();System.out.println(pair.atom1().toString()+ "  " + pair.atom2().toString());}  line = in.readLine();
 
        System.out.println("reset(atom,BOTH)"); iterator.reset(id.set(atom).set(IteratorDirective.BOTH));
        while(iterator.hasNext()) {AtomPair pair = iterator.next();System.out.println(pair.atom1().toString()+ "  " + pair.atom2().toString());}  line = in.readLine();
 
        System.out.println("reset(UP)"); iterator.reset(id.set().set(IteratorDirective.UP));
        while(iterator.hasNext()) {AtomPair pair = iterator.next();System.out.println(pair.atom1().toString()+ "  " + pair.atom2().toString());}  line = in.readLine();
 
        System.out.println("reset(DOWN)"); iterator.reset(id.set().set(IteratorDirective.DOWN));
        while(iterator.hasNext()) {AtomPair pair = iterator.next();System.out.println(pair.atom1().toString()+ "  " + pair.atom2().toString());}  line = in.readLine();
 
        System.out.println("reset(NEITHER)"); iterator.reset(id.set().set(IteratorDirective.NEITHER));
        while(iterator.hasNext()) {AtomPair pair = iterator.next();System.out.println(pair.atom1().toString()+ "  " + pair.atom2().toString());}  line = in.readLine();
 
        System.out.println("reset(BOTH)"); iterator.reset(id.set().set(IteratorDirective.BOTH));
        while(iterator.hasNext()) {AtomPair pair = iterator.next();System.out.println(pair.atom1().toString()+ "  " + pair.atom2().toString());}  line = in.readLine();
        System.out.println("Done");
        line = in.readLine();
        System.exit(0);
    }//end of main
 */   
}  //end of class AtomPairIterator
    
