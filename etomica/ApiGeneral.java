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
     * The iterators used to generate the sets of atoms.
     * The inner one is not necessarily atom dependent.
     */
    private AtomIterator aiInner, aiOuter;
    
    protected boolean hasNext;

    
    /**
     * Construct a pair iterator using the given atom iterators
     */
    public ApiGeneral(Space s, AtomIterator aiOuter, AtomIterator aiInner) {
        pair = new AtomPair(s);
        hasNext = false;
        this.aiOuter = aiOuter;
        this.aiInner = aiInner;
    }
    
    public AtomIterator aiOuter() {return aiOuter;}
    public AtomIterator aiInner() {return aiInner;}

    public void setInnerIterator(AtomIterator inner) {
    	this.aiInner = inner;
    	unset();
    }
    
    public void setOuterIterator(AtomIterator outer) {
    	this.aiOuter = outer;
    	unset();
    }
    
    public void unset() {
    	hasNext = false;
    }
    
    /**
     * Returns whether or not an AtomPair is contained by the
     * iterator. The first atom of the pair must be contained in
     * the outer loop.  The second atom of the pair must be contained
     * in the inner loop.
     */
    public boolean contains(AtomPair pair) {
    	return aiOuter.contains(pair.atom1) && aiInner.contains(pair.atom2());
    }
    
    
    /**
     * Returns the number of pairs given by this iterator.
     */
    public int size() {
    	return aiOuter.size() * aiInner.size();
    }        
    
    public final boolean hasNext() {return hasNext;}

    /**
     * Resets the iterator, so that it is ready to go through all of its pairs.
     */
    public void reset() {
        aiOuter.reset();
        aiInner.reset();
        hasNext = aiOuter.hasNext() && aiInner.hasNext();
    }
    
    /**
     * Returns the next pair without advancing the iterator.
     * If the iterator has reached the end of its iteration,
     * returns null.
     */
    public AtomPair peek() {
    	if(!hasNext) {return null;}
    	
    	if(aiInner.hasNext()) {
    		pair.reset2(aiInner.peek());
    	}
    	else {
    		// We can reset the inner loop, because it has
    		// reached its end, and the next() method
    		// will see the inner loop at its
    		// start, and advance from there.
    		aiInner.reset();
    		pair.reset(aiOuter.next(), aiInner.peek());
    	}
    	return pair;
    }
    
    /**
     * Returns the next pair of atoms.  Advances the outer loop to 
     * the next atom when the inner loop has reached its last atom.
     */
    public final AtomPair next() {
    	if(!hasNext) {return null;}
    	//Advance the inner loop, if it is not at its end.
    	if(aiInner.hasNext()) {
    		pair.reset2(aiInner.next());
    	}
    	//Advance the outer loop, if the inner loop has reached its end.
    	else {
    		aiInner.reset();
			pair.reset(aiOuter.next() , aiInner.next());
    	} 
    	
    	hasNext = aiInner.hasNext() || aiOuter.hasNext();
        return pair;
    }

    /**
     * Performs the given action on all pairs returned by this iterator.
     */
    public void allPairs(AtomPairActive act) {  
       reset();
       while(aiOuter.hasNext()) {
       		pair.setAtom1(aiOuter.next());
       		aiInner.reset();
       		while(aiInner.hasNext()){
       			act.actionPerformed(pair.reset2(aiInner.next()));
       		}
       }
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
    
