package etomica;

/**
 * Basic class for iterating over triples of atoms.
 * Triplets are iterated by collecting atoms yielded by three atom iterators.
 * Different types of triplet iterators can be constructed with different choices
 * of the atom iterators.
 *
 * @author David Kofke
 */
public class Atom3Iterator implements java.io.Serializable {
    
    private Atom3 atom3; //want final, but Null inner class won't allow
    private IteratorDirective.Direction direction;
    private final IteratorDirective localDirective = new IteratorDirective();

    /**
     * The iterators used to generate the sets of atoms
     */
    protected /*final*/ AtomIterator ai1, ai2, ai3;
    protected AtomIterator aiInner, aiOuter, aiMiddle;
    
    /**
     * A pair action wrapper used to enable the allPairs method
     */
 //   protected AtomPairAction.Wrapper actionWrapper;   // want final too //wrapper inner class defined below
    protected boolean hasNext;
    /**
     * Flag indicating whether atom1 of pair needs to be updated to point to the same atom that "atom1" in this class points to
     */
    private boolean needUpdate1; 
    private Atom atom1;
    
    //Used only for the NULL iterator, defined below.
    private Atom3Iterator() {
        hasNext = false;
        pair = null;
        ai1 = ai2 = ai3 = AtomIterator.NULL;
    }
    /**
     * Constructs an iterator of all atom pairs in the given phase.
     */
    public Atom3Iterator(Phase p) {
        this(p.parentSimulation().space, p.iteratorFactory().makeAtomIterator(), p.iteratorFactory().makeAtomIterator());
    }
    /**
     * Constructs an iterator of all pairs formed from the given species in the given phase.
     */
/*     public AtomPairIterator(Phase p, Species species1, Species species2) {
        this(p, species1.getAgent(p).new LeafAtomIterator(),
                species2.getAgent(p).new LeafAtomIterator());
     }*/
     /** 
      * Sets pair iterator so that it traverses all leaf-atom pairs in its basis.
      */
     public Atom3Iterator(Space s) {
        this(s, new AtomIteratorSequential(true), new AtomIteratorSequential(true));
     }
    /**
     * Construct a triplet iterator using the given atom iterators
     */
    public Atom3Iterator(Space s, AtomIterator iter1, AtomIterator iter2, AtomItertor iter3) {
        atom3 = new Atom3(s);
        hasNext = false;
        ai1 = iter1;
        ai2 = iter2;
        ai3 = iter3;
        direction = IteratorDirective.UP;
    }
    
    public void setBasis(Atom a1, Atom a2, Atom3 a3) {
        ai1.setBasis(a1);
        ai2.setBasis(a2);
        ai3.setBasis(a3);
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
    
    private final void setOuterInner(AtomIterator outIter, AtomIterator midIter, AtomIterator inIter) {
        aiOuter = outIter;
        aiMiddle = midIter;
        aiInner = inIter;
        aiOuter.setAsNeighbor(false);
        aiMiddle.setAsNeighbor(true);
        aiInner.setAsNeighbor(true);
    }

    /**
     * Resets the iterator, so that it is ready to go through all of its triplets.
     */
    public void reset() {
        if(direction == IteratorDirective.BOTH) direction = IteratorDirective.UP;
        if(ai2.getBasis().preceeds(ai1.getBasis())) {
            if(direction.doUp()) setOuterInner(ai2, ai1);
            else setOuterInner(ai1, ai2);
        } else {
            if(direction.doUp()) setOuterInner(ai1, ai2);
            else setOuterInner(ai2, ai1);
        }
    //    setOuterInner(ai1, ai2);
        aiOuter.reset(localDirective.set().set(direction));
        aiInner.reset(localDirective);
        setFirst();
    }
        
    /**
     * Resets the iterator so that it iterates over all pairs formed with the 
     * given atom in the most recently specified iterator directive (default UP is
     * if none previously specified.
     */
    public void reset(Atom atom) {
        if(ai1.contains(atom)) setOuterInner(ai1, ai2);
        else if(ai2.contains(atom)) setOuterInner(ai2, ai1);
        else {hasNext = false; return;}
        aiOuter.reset(localDirective.set(atom).set(IteratorDirective.NEITHER));
        aiInner.reset(localDirective.set().set(direction));//reset inner last to leave localDirective in proper state for inner-loop resets
        setFirst();
    }
    
    /**
     * This is called after aiOuter, aiMiddle and aiInner are are reset, 
     * and it readies iterators to give first pair (or puts hasNext = false if no pairs
     * are forthcoming).
     */
    protected final void setFirst() {
        hasNext = false;
        while(aiOuter.hasNext()) { //loop over iterator 1...
            pair.atom1 = aiOuter.next();
            aiInner.reset(localDirective.set(pair.atom1));
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
        if(needUpdate2) {pair.atom2 = aimt2; needUpdate2 = false;}  //aiMiddle was advanced
        atom3.atom3 = aiInner.next();
        atom3.reset();
        while(!aiMiddle.hasNext()) {//Middle is done for this atom1, loop until it is prepared for next
            if(aiOuter.hasNext()) {     //Outer has another atom1...
                atom1 = aiOuter.next();
                aiMiddle.reset(localDirective.set(atom1));
                needUpdate1 = true;
                
            
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
     * Performs the given action on all pairs returned by this iterator.
     */
     //not carefully checked
/*    public void allPairs(AtomPairAction act) {  
        reset();
//        ai1.reset();  //this shouldn't be here, in general; need to consider it more carefully
        actionWrapper.pairAction = act;
        while(ai1.hasNext()) {
            pair.atom1 = ai1.next();
            ai2.reset(localDirective.set(pair.atom1));
            ai2.allAtoms(actionWrapper);
        }
    }
*/    
    
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
    
