package etomica;

/**
 * Encapsulation of a set of instructions that an (Atom/AtomPair)Iterator
 * uses to select the atoms/pairs it presents on iteration.
 *
 * @author David Kofke
 */
public class IteratorDirective implements java.io.Serializable {
    
    private Atom atom1, atom2;
    private Direction direction;
    private int atomCount;
    private PotentialCriterion potentialCriteriaHead;
    
    public IteratorDirective() {
        this(UP);
    }
    public IteratorDirective(Direction direction) {
        set(direction);
        set();
    }
    
    //returns itself as a convenience, so that it may be set while being passed as an
    //argument to a method
    public final IteratorDirective set() {
        atom1 = atom2 = null;
        atomCount = 0;
        return this;
    }
    public final IteratorDirective set(Atom a) {
        atom1 = a;
        atom2 = null;
        atomCount = 1;
        return this;
    }
    public final IteratorDirective set(Atom a1, Atom a2) {
        atom1 = a1;
        atom2 = a2;
        atomCount = 2;
        return this;
    }
    public final IteratorDirective set(Direction direction) {
        this.direction = direction;
        return this;
    }

    public final int atomCount() {return atomCount;}
    public final Direction direction() {return direction;}
    
    public final Atom atom1() {return atom1;}
    public final Atom atom2() {return atom2;}
    
    public final boolean excludes(PotentialAgent p) {
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
     */
    public static abstract class PotentialCriterion {
        /**
         * Definition of criterion.  If this method returns true, the potential's atoms
         * are excluded from iteration.
         */
        public abstract boolean excludes(PotentialAgent pot);
        
        //Linked-list constructs
        private PotentialCriterion nextCriterion;
        public void setNextCriterion(PotentialCriterion next) {nextCriterion = next;}
        public PotentialCriterion nextCriterion() {return nextCriterion;}
    }//end of PotentialCriterion
    
}//end of IteratorDirective    