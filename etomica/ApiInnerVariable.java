package etomica;

import etomica.action.AtomsetAction;

/**
 * Pair iterator synthesized from two atom iterators, such that the inner-loop
 * iteration depends on the outer-loop atom.  Pairs are formed from
 * the atoms yielded by the two atom iterators.  The inner-loop iterator must
 * implement AtomIteratorAtomDependent, and its set(Atom) method will be invoked
 * with the current outer-loop atom at the start of each inner-loop iteration.
 * All pairs returned by iterator are the same AtomPair instance, and
 * differ only in the Atom instances held by it.
 */
 
 /* History of changes
  * 08/25/04 (DAK et al) new
  */
  
public class ApiInnerVariable implements AtomsetIterator, ApiComposite {
    
    /**
     * Construct a pair iterator using the given atom iterators.  Requires
     * call to reset() before beginning iteration.
     */
    public ApiInnerVariable(AtomIterator aiOuter, 
    						AtomIteratorAtomDependent aiInner) {
        this.aiOuter = aiOuter;
        this.aiInner = aiInner;
        unset();
    }
    
    /**
     * Accessor method for the outer-loop atom iterator.
     * @return the current outer-loop iterator
     */
    public AtomIterator getOuterIterator() {return aiOuter;}
    
    /**
     * Accessor method for the inner-loop atom iterator.
     * @return the current inner-loop iterator
     */
    public AtomIterator getInnerIterator() {return aiInner;}

    /**
     * Defines the atom iterator that performs the inner-loop iteration to
     * generate the pairs.
     * @param inner The new inner-loop iterator.
     */
    public void setInnerIterator(AtomIteratorAtomDependent inner) {
    	this.aiInner = inner;
    	unset();
    }
    
    /**
     * Defines the iterator the performs the outer-loop iteration to
     * generate the pairs.
     * @param outer The new outer-loop iterator.
     */
    public void setOuterIterator(AtomIterator outer) {
    	this.aiOuter = outer;
    	unset();
    }
    
    /**
     * Sets the iterator such that hasNext is false.
     */
    public void unset() {
    	hasNext = false;
    }
    
    /**
     * Indicates whether the given atom pair will be returned by the
     * iterator during its iteration. The order of the atoms in the pair
     * is significant (this means that a value of true is returned only if
     * one of the pairs returned by the iterator will have the same two 
     * atoms in the same atom1/atom2 position as the given pair). Not dependent
     * on state of hasNext.
     */
    public boolean contains(Atom[] pair) {
    	if(aiOuter.contains(new Atom[] {pair[0]})) {
    		aiInner.setAtom(pair[0]);
    		return aiInner.contains(new Atom[] {pair[1]});
    	}
        return false;
    }

    /**
     * Returns the number of pairs given by this iterator.  Independent
     * of state of hasNext. Clobbers the iteration state (i.e., status
     * of hasNext/next) but does not recondition iterator (i.e., does not 
     * change set of iterates that would be given on iteration after reset).
     * Must perform reset if attempting iteration after using size() method.
     */
    public int size() {
    	int sum = 0;
        aiOuter.reset();
        while(aiOuter.hasNext()) {
    		aiInner.setAtom(aiOuter.nextAtom());
    		sum += aiInner.size();
        }
    	return sum;
    }        
    
    /**
     * Indicates whether the iterator has completed its iteration.
     */
    public boolean hasNext() {return hasNext;}

    /**
     * Resets the iterator, so that it is ready to go through all of its pairs.
     */
    public void reset() {
        aiOuter.reset();
        hasNext = false;
        needUpdate1 = false;
        while(aiOuter.hasNext()) { //loop over outer iterator...
            pair[0] = aiOuter.nextAtom();
            aiInner.setAtom(pair[0]);
            aiInner.reset();
            if(aiInner.hasNext()) { //until inner iterator has another
                hasNext = true;
                break;        //...until iterator 2 hasNext
            }
        }//end while
    }
    
    /**
     * Returns the next pair without advancing the iterator.
     * If the iterator has reached the end of its iteration,
     * returns null.
     */
    public Atom[] peek() {
    	if(!hasNext) {return null;}   	
        if(needUpdate1) {pair[0] = atom1; needUpdate1 = false;}  //aiOuter was advanced
        pair[1] = aiInner.peek()[0];
    	return pair;
    }
    
    /**
     * Returns the next pair of atoms. The same AtomPair instance
     * is returned every time, but the Atoms it holds are (of course)
     * different for each iterate. 
     */
    public Atom[] next() {
    	if(!hasNext) return null;
        //we use this update flag to indicate that atom1 in pair needs to be set to a new value.
        //it is not done directly in the while-loop because pair must first return with the old atom1 intact
        if(needUpdate1) {pair[0] = atom1; needUpdate1 = false;}  //aiOuter was advanced
        pair[1] = aiInner.nextAtom();
        while(!aiInner.hasNext()) {     //Inner is done for this atom1, loop until it is prepared for next
            if(aiOuter.hasNext()) {     //Outer has another atom1...
                atom1 = aiOuter.nextAtom();           //...get it
                aiInner.setAtom(atom1);
                aiInner.reset();
                needUpdate1 = true;           //...flag update of pair.atom1 for next time
            }
            else {hasNext = false; break;} //Outer has no more; all done with pairs
        }//end while
        return pair;
    }

    /**
     * Performs the given action on all pairs returned by this iterator.
     */
    public void allAtoms(AtomsetAction act) {  
		aiOuter.reset();
		while(aiOuter.hasNext()) {
			pair[0] = aiOuter.nextAtom();
			aiInner.setAtom(pair[0]);
			aiInner.reset();
			while(aiInner.hasNext()){
				pair[1] = aiInner.nextAtom();
				act.actionPerformed(pair);
			}
		}
    }
    
    public final int nBody() {return 2;}
    
    protected final Atom[] pair = new Atom[2];
    protected boolean hasNext, needUpdate1;
    protected Atom atom1;

    /**
     * The iterators used to generate the sets of atoms.
     * The inner one is not necessarily atom dependent.
     */
    protected AtomIteratorAtomDependent aiInner;
	protected AtomIterator aiOuter;
    
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
    
