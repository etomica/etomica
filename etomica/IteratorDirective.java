package etomica;

/**
 * Encapsulation of a set of instructions that an (Atom/AtomPair)Iterator
 * uses to select the atoms/pairs it presents on iteration.
 *
 * @author David Kofke
 */
public final class IteratorDirective implements java.io.Serializable {
    
    private Atom atom1;
    private Direction direction;
    private int atomCount;
    PotentialCriterion potentialCriteriaHead;
    
    /**
     * Flag indicating whether long-range correction contributions should
     * be included in calculation.  Default is <b>false</b>.
     */
    public boolean includeLrc = false;
    
    public IteratorDirective() {
        this(UP);
    }
    public IteratorDirective(Direction direction) {
        set(direction);
        set();
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
     */
    public void copy(IteratorDirective id) {
        direction = id.direction();
        includeLrc = id.includeLrc;
        atom1 = id.atom1();
        atomCount = (atom1 == null) ? 0 : 1;
        for(PotentialCriterion crit=id.potentialCriteriaHead; crit!=null; crit=crit.nextCriterion()) {
            addCriterion((PotentialCriterion)crit.clone());
        }
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
    
    /**
     * Method for testing of iterators.  The iterates of the given iterator are generated
     * and their signature is written to the console.  Various combinations of iteration direction
     * and reference atom are tested.  The given atoms are used for the reference in
     * some reset calls, and permit examintion of the iterator's behavior when the first and
     * last atoms are used for the reference, and an atom chosen from the middle of its
     * sequence.  Both hasNext/next and allAtoms methods of iteration are performed.
     */
    public static void testSuite(AtomIterator iterator, Atom first, Atom middle, Atom last) {
        
        IteratorDirective directive = new IteratorDirective();
        try {
            System.out.println("reset()");
            iterator.reset();
            test(iterator);
        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); /*pauseForInput();*/}
        
        try {
            System.out.println("reset(directive.set(middle).set(IteratorDirective.BOTH))");
            iterator.reset(directive.set(middle).set(IteratorDirective.BOTH));
            test(iterator);
        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); /*pauseForInput();*/}
            
        try {
            System.out.println("reset(directive.set(middle).set(IteratorDirective.UP))");
            iterator.reset(directive.set(middle).set(IteratorDirective.UP));
            test(iterator);
        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
            
        try {
            System.out.println("reset(directive.set(middle).set(IteratorDirective.DOWN))");
            iterator.reset(directive.set(middle).set(IteratorDirective.DOWN));
            test(iterator);
        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
            
        try {
            System.out.println("reset(directive.set(first).set(IteratorDirective.BOTH))");
            iterator.reset(directive.set(first).set(IteratorDirective.BOTH));
            test(iterator);
        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
        try {
            System.out.println("reset(directive.set(first).set(IteratorDirective.UP))");
            iterator.reset(directive.set(first).set(IteratorDirective.UP));
            test(iterator);
        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
        try {
            System.out.println("reset(directive.set(first).set(IteratorDirective.DOWN))");
            iterator.reset(directive.set(first).set(IteratorDirective.DOWN));
            test(iterator);
        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
            
        try {
            System.out.println("reset(directive.set(last).set(IteratorDirective.BOTH))");
            iterator.reset(directive.set(last).set(IteratorDirective.BOTH));
            test(iterator);
        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
        try {
            System.out.println("reset(directive.set(last).set(IteratorDirective.UP))");
            iterator.reset(directive.set(last).set(IteratorDirective.UP));
            test(iterator);
        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
        try {
            System.out.println("reset(directive.set(last).set(IteratorDirective.DOWN))");
            iterator.reset(directive.set(last).set(IteratorDirective.DOWN));
            test(iterator);
        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
            
        try {
            System.out.println("reset(directive.set().set(IteratorDirective.BOTH))");
            iterator.reset(directive.set().set(IteratorDirective.BOTH));
            test(iterator);
        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
        try {
            System.out.println("reset(directive.set().set(IteratorDirective.UP))");
            iterator.reset(directive.set().set(IteratorDirective.UP));
            test(iterator);
        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
        try {
            System.out.println("reset(directive.set().set(IteratorDirective.DOWN))");
            iterator.reset(directive.set().set(IteratorDirective.DOWN));
            test(iterator);
        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
            
        try {
            System.out.println("reset(directive.set(first).set(IteratorDirective.NEITHER))");
            iterator.reset(directive.set(first).set(IteratorDirective.NEITHER));
            test(iterator);
        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
            
        try {
            System.out.println("reset(directive.set().set(IteratorDirective.NEITHER))");
            iterator.reset(directive.set().set(IteratorDirective.NEITHER));
            test(iterator);
        } catch(RuntimeException ex) {System.out.println("RuntimeException caught: "+ex.getMessage()); pauseForInput();}
    }
        
    /**
     * Causes signature of all iterates of the given iterator, as currently reset,
     * to be printed to console.  Performs iteration using hasNext/next loop, then
     * using allAtoms method (done twice).  Support method for testSuite method.
     */
    private static void test(AtomIterator iterator) {
        AtomAction printAtom = new AtomAction() {
            public void actionPerformed(Atom a) {System.out.println(a.signature());}
        };
        
        while(iterator.hasNext()) System.out.println(iterator.next().signature());
        pauseForInput();
        System.out.println("allAtoms (twice)");
        iterator.allAtoms(printAtom);
        System.out.println();
        iterator.allAtoms(printAtom);
        pauseForInput();
    }
    /**
     * Halts program activity until a return is entered from the console.
     * Support method for testSuite method.
     */
    private static void pauseForInput() {
        System.out.println("Hit return to continue");
        try {
            System.in.read();
            System.in.read();
        } catch(Exception e) {}
    }
    
    //IteratorDirective.Direction
    public static final class Direction extends Constants.TypedConstant {
        
        private final boolean up;
        private final boolean down;
        private Direction(String label, boolean doUp, boolean doDown) {
            super(label);
            up = doUp;
            down = doDown;
        }
        public boolean doUp() {return up;}
        public boolean doDown() {return down;}
        
        public static final Direction[] CHOICES = new Direction[] {
            new Direction("Up", true, false),
            new Direction("Down", false, true),
            new Direction("Neither", false, false),
            new Direction("Both", true, true)
        };
        
        public final Constants.TypedConstant[] choices() {return CHOICES;}
    }//end of Direction
    public static final Direction UP = Direction.CHOICES[0];
    public static final Direction DOWN = Direction.CHOICES[1];
    public static final Direction NEITHER = Direction.CHOICES[2];
    public static final Direction BOTH = Direction.CHOICES[3];
    
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
    
}//end of IteratorDirective    