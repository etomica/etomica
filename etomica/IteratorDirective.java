package etomica;

/**
 * Encapsulation of a set of instructions that an AtomsetIterator
 * uses to select the atoms it presents on iteration.
 *
 * @author David Kofke
 */

/* History
 * 02/21/03 (DAK) added constructor taking Direction and Atom
 * 08/30/03 (DAK) added copying of skipFirst to copy method
 * 
 */
public final class IteratorDirective implements java.io.Serializable {
    
    private Atom atom1;
    private Direction direction;
    private int atomCount;
    PotentialCriterion potentialCriteriaHead;
    
    /**
     * Flag indicating whether long-range correction contributions should
     * be included in calculation.  Default is <b>true</b>.
     */
    public boolean includeLrc = true;
        
    public IteratorDirective() {
        this(UP);
    }
    public IteratorDirective(Direction direction) {
        set(direction);
        set();
    }
    public IteratorDirective(Direction direction, Atom atom) {
    	set(direction);
    	set(atom);
    }
    
    /**
     * Puts directive in default state of no atoms specified, up direction, no
     * potential criteria applied, no LRC included.
     */
    public IteratorDirective clear() {
        set(UP);
        set();
        potentialCriteriaHead = null;
        includeLrc = false;
        return this;
    }
    
    /**
     * Puts all settings of this directive to equal those of the given directive.
     * The list of potential criteria of the given directive is referenced as the list
     * for this directive, so changes in list for either directive will affect the
     * list for the other directive.
     */
     //we don't make copy of list to avoid overhead of construction
    public void copy(IteratorDirective id) {
        direction = id.direction();
        includeLrc = id.includeLrc;
        atom1 = id.atom1();
        atomCount = (atom1 == null) ? 0 : 1;
        potentialCriteriaHead = id.potentialCriteriaHead;
    }
    
    //returns itself as a convenience, so that it may be set while being passed as an
    //argument to a method
    public final IteratorDirective set() {
        atom1 = null;
        atomCount = 0;
        return this;
    }
    public final IteratorDirective set(Atom a) {
        atom1 = a;
        atomCount = (atom1 != null) ? 1 : -1;
        return this;
    }
    public final IteratorDirective set(Direction direction) {
        this.direction = direction;
        return this;
    }

    public final int atomCount() {return atomCount;}
    public final Direction direction() {return direction;}
    
    public final Atom atom1() {return atom1;}
    
    public final boolean excludes(Potential p) {
        for(PotentialCriterion crit=potentialCriteriaHead; crit!=null; crit=crit.nextCriterion()) {
            if(crit.excludes(p)) return true;
        }
        return false;
    }
    
    /**
     * Adds criterion to set of criteria for potentials.  Criteria
     * are considered in the reverse order of their addition (last-added is considered first).
     * There is no way to remove a criterion.
     */
    public final void addCriterion(PotentialCriterion newCriterion) {
        if(newCriterion == null) return;
        newCriterion.setNextCriterion(potentialCriteriaHead);
        potentialCriteriaHead = newCriterion;
    }
    
    public void setTargetAtoms(Atom[] atoms) {targetAtoms = atoms;}
    public Atom[] targetAtoms() {return targetAtoms;}
    
    public Atom[] targetAtoms;
    
    /**
     * Method for testing of iterators.  The iterates of the given iterator are generated
     * and their signature is written to the console.  Various combinations of iteration direction
     * and reference atom are tested.  The given atoms are used for the reference in
     * some reset calls, and permit examination of the iterator's behavior when the first and
     * last atoms are used for the reference, and an atom chosen from the middle of its
     * sequence.  Both hasNext/next and allAtoms methods of iteration are performed.
     */
//    public static void testSuite(AtomIterator iterator, Atom first, Atom middle, Atom last) {
//        
//        IteratorDirective directive = new IteratorDirective();
//        try {
//            System.out.println("reset()");
//            iterator.reset();
//            test(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); /*pauseForInput();*/}
//        
//        try {
//            System.out.println("reset(directive.set(middle).set(IteratorDirective.BOTH))");
//            iterator.reset(directive.set(middle).set(IteratorDirective.BOTH));
//            test(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); /*pauseForInput();*/}
//            
//        try {
//            System.out.println("reset(directive.set(middle).set(IteratorDirective.UP))");
//            iterator.reset(directive.set(middle).set(IteratorDirective.UP));
//            test(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
//            
//        try {
//            System.out.println("reset(directive.set(middle).set(IteratorDirective.DOWN))");
//            iterator.reset(directive.set(middle).set(IteratorDirective.DOWN));
//            test(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
//            
//        try {
//            System.out.println("reset(directive.set(first).set(IteratorDirective.BOTH))");
//            iterator.reset(directive.set(first).set(IteratorDirective.BOTH));
//            test(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
//        try {
//            System.out.println("reset(directive.set(first).set(IteratorDirective.UP))");
//            iterator.reset(directive.set(first).set(IteratorDirective.UP));
//            test(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
//        try {
//            System.out.println("reset(directive.set(first).set(IteratorDirective.DOWN))");
//            iterator.reset(directive.set(first).set(IteratorDirective.DOWN));
//            test(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
//            
//        try {
//            System.out.println("reset(directive.set(last).set(IteratorDirective.BOTH))");
//            iterator.reset(directive.set(last).set(IteratorDirective.BOTH));
//            test(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
//        try {
//            System.out.println("reset(directive.set(last).set(IteratorDirective.UP))");
//            iterator.reset(directive.set(last).set(IteratorDirective.UP));
//            test(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
//        try {
//            System.out.println("reset(directive.set(last).set(IteratorDirective.DOWN))");
//            iterator.reset(directive.set(last).set(IteratorDirective.DOWN));
//            test(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
//            
//        try {
//            System.out.println("reset(directive.set().set(IteratorDirective.BOTH))");
//            iterator.reset(directive.set().set(IteratorDirective.BOTH));
//            test(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
//        try {
//            System.out.println("reset(directive.set().set(IteratorDirective.UP))");
//            iterator.reset(directive.set().set(IteratorDirective.UP));
//            test(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
//        try {
//            System.out.println("reset(directive.set().set(IteratorDirective.DOWN))");
//            iterator.reset(directive.set().set(IteratorDirective.DOWN));
//            test(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
//            
//        try {
//            System.out.println("reset(directive.set(first).set(IteratorDirective.NEITHER))");
//            iterator.reset(directive.set(first).set(IteratorDirective.NEITHER));
//            test(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
//            
//        try {
//            System.out.println("reset(directive.set().set(IteratorDirective.NEITHER))");
//            iterator.reset(directive.set().set(IteratorDirective.NEITHER));
//            test(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
//    }
//        
//    /**
//     * Causes signature of all iterates of the given iterator, as currently reset,
//     * to be printed to console.  Performs iteration using hasNext/next loop, then
//     * using allAtoms method (done twice).  Support method for testSuite method.
//     */
//    private static void test(AtomIterator iterator) {
//        AtomAction printAtom = new AtomAction() {
//            public void actionPerformed(Atom a) {System.out.println(a.signature());}
//        };
//        
//        while(iterator.hasNext()) System.out.println(iterator.next().signature());
//        pauseForInput();
//        System.out.println("allAtoms (twice)");
//        iterator.allAtoms(printAtom);
//        System.out.println();
//        iterator.allAtoms(printAtom);
//        pauseForInput();
//    }
//    
//    /**
//     * Method for testing of pair iterators.  The iterates of the given iterator are generated
//     * and their signature is written to the console.  Various combinations of iteration direction
//     * and reference atom are tested.  The given atoms are used for the reference in
//     * some reset calls, and permit examination of the iterator's behavior when the first and
//     * last atoms are used for the reference, and an atom chosen from the middle of its
//     * sequence.  Both hasNext/next and allAtoms methods of iteration are performed.
//     */
//    public static void testSuitePair(AtomPairIterator iterator, Atom first, Atom middle, Atom last) {
//        
//        IteratorDirective directive = new IteratorDirective();
//        try {
//            System.out.println("reset()");
//            iterator.reset();
//            testPair(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); /*pauseForInput();*/}
//        
//        try {
//            System.out.println("reset(directive.set(middle).set(IteratorDirective.BOTH))");
//            iterator.reset(directive.set(middle).set(IteratorDirective.BOTH));
//            testPair(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); /*pauseForInput();*/}
//            
//        try {
//            System.out.println("reset(directive.set(middle).set(IteratorDirective.UP))");
//            iterator.reset(directive.set(middle).set(IteratorDirective.UP));
//            testPair(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
//            
//        try {
//            System.out.println("reset(directive.set(middle).set(IteratorDirective.DOWN))");
//            iterator.reset(directive.set(middle).set(IteratorDirective.DOWN));
//            testPair(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
//            
//        try {
//            System.out.println("reset(directive.set(first).set(IteratorDirective.BOTH))");
//            iterator.reset(directive.set(first).set(IteratorDirective.BOTH));
//            testPair(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
//        try {
//            System.out.println("reset(directive.set(first).set(IteratorDirective.UP))");
//            iterator.reset(directive.set(first).set(IteratorDirective.UP));
//            testPair(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
//        try {
//            System.out.println("reset(directive.set(first).set(IteratorDirective.DOWN))");
//            iterator.reset(directive.set(first).set(IteratorDirective.DOWN));
//            testPair(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
//            
//        try {
//            System.out.println("reset(directive.set(last).set(IteratorDirective.BOTH))");
//            iterator.reset(directive.set(last).set(IteratorDirective.BOTH));
//            testPair(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
//        try {
//            System.out.println("reset(directive.set(last).set(IteratorDirective.UP))");
//            iterator.reset(directive.set(last).set(IteratorDirective.UP));
//            testPair(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
//        try {
//            System.out.println("reset(directive.set(last).set(IteratorDirective.DOWN))");
//            iterator.reset(directive.set(last).set(IteratorDirective.DOWN));
//            testPair(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
//            
//        try {
//            System.out.println("reset(directive.set().set(IteratorDirective.BOTH))");
//            iterator.reset(directive.set().set(IteratorDirective.BOTH));
//            testPair(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
//        try {
//            System.out.println("reset(directive.set().set(IteratorDirective.UP))");
//            iterator.reset(directive.set().set(IteratorDirective.UP));
//            testPair(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
//        try {
//            System.out.println("reset(directive.set().set(IteratorDirective.DOWN))");
//            iterator.reset(directive.set().set(IteratorDirective.DOWN));
//            testPair(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
//            
//        try {
//            System.out.println("reset(directive.set(first).set(IteratorDirective.NEITHER))");
//            iterator.reset(directive.set(first).set(IteratorDirective.NEITHER));
//            testPair(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
//            
//        try {
//            System.out.println("reset(directive.set().set(IteratorDirective.NEITHER))");
//            iterator.reset(directive.set().set(IteratorDirective.NEITHER));
//            testPair(iterator);
//        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
//        
//    }
//    
//    /**
//     * Causes signature of all pair iterates of the given iterator, as currently reset,
//     * to be printed to console.  Performs iteration using hasNext/next loop, then
//     * using allAtoms method (done twice).  Support method for testSuitePair method.
//     */
//    private static void testPair(AtomPairIterator iterator) {
//        AtomPairAction printAtom = new AtomPairAction(Simulation.instance.space) {
//            public void actionPerformed(AtomPair pair) {System.out.println(pair.atom1().signature()+" "+pair.atom2().signature());}
//        };
//        
//        while(iterator.hasNext()) {
//            AtomPair pair = iterator.next();
//            System.out.println(pair.atom1().signature()+" "+pair.atom2().signature());
//        }
//        pauseForInput();
//        System.out.println("allAtoms (twice)");
//        iterator.allPairs(printAtom);
//        System.out.println();
//        iterator.allPairs(printAtom);
//
//        pauseForInput();
//    }
//    
//    
    /**
     * Halts program activity until a return is entered from the console.
     * Support method for testSuite method.
     */
    public static void pauseForInput() {
        System.out.println("Hit return to continue");
        try {
            System.in.read();
            System.in.read();
        } catch(Exception e) {}
    }
    
    //IteratorDirective.Direction
    public static final class Direction extends Constants.TypedConstant {
        
        private Direction(String label) {
            super(label);
        }
        public static final Direction[] CHOICES = new Direction[] {
            new Direction("Up"),
            new Direction("Down"),
        };
        
        public final Constants.TypedConstant[] choices() {return CHOICES;}
    }//end of Direction
    public static final Direction UP = Direction.CHOICES[0];
    public static final Direction DOWN = Direction.CHOICES[1];
    
    /**
     * Class used to define a criterion that must be satisfied by a potential
     * in order for its atoms to be approved for iteration by an iterator directive.
     * Multiple criteria are ordered into a linked list by the iterator directive.
     * This is made cloneable to support IteratorDirective.copy functionality.
     */
    public static abstract class PotentialCriterion implements Cloneable {
        /**
         * Definition of criterion.  If this method returns true, the potential's atoms
         * are excluded from iteration.
         * Note that any subclasses should be sure to override clone method if more than
         * a shallow copy is appropriate.
         */
        public abstract boolean excludes(Potential pot);
        
        //Linked-list constructs
        private PotentialCriterion nextCriterion;
        public void setNextCriterion(PotentialCriterion next) {nextCriterion = next;}
        public PotentialCriterion nextCriterion() {return nextCriterion;}
        
        public Object clone() {
            Object obj = null;
            try {
                obj = super.clone();
            } catch(CloneNotSupportedException e) {e.printStackTrace();}
            return obj;
        }
    }//end of PotentialCriterion
    
	/**
	 * Sets flag indicating if lrc potentials (long-range correction) should be
	 * included.
	 * @return boolean
	 */
	public boolean isIncludeLrc() {
		return includeLrc;
	}

	/**
	 * Sets flag indicating if lrc potentials (long-range correction) should be
	 * included.
	 * @param includeLrc The flag value to set
	 * @return this IteratorDirective, for in-line use of the method.
	 */
	public IteratorDirective setIncludeLrc(boolean includeLrc) {
		this.includeLrc = includeLrc;
		return this;
	}


}//end of IteratorDirective    