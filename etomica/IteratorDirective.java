package etomica;

/**
 * Encapsulation of a set of instructions that an (Atom/AtomPair)Iterator
 * uses to select the atoms/pairs it presents on iteration.
 *
 * @author David Kofke
 */
public final class IteratorDirective implements java.io.Serializable {
    
    private Atom atom1, atom2;
    private Direction direction;
    private int atomCount;
    PotentialCriterion potentialCriteriaHead;
    
    public IteratorDirective() {
        this(UP);
    }
    public IteratorDirective(Direction direction) {
        set(direction);
        set();
    }
    
    /**
     * Puts all settings of this directive to equal those of the given directive.
     */
    public void copy(IteratorDirective id) {
        direction = id.direction();
        atom1 = id.atom1();
        atom2 = id.atom2();
        if(atom1 == null) atomCount = 0;
        else if(atom2 == null) atomCount = 1;
        else atomCount = 2;
        for(PotentialCriterion crit=id.potentialCriteriaHead; crit!=null; crit=crit.nextCriterion()) {
            addCriterion((PotentialCriterion)crit.clone());
        }
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
        atomCount = (atom1 != null) ? 1 : -1;
        return this;
    }
    public final IteratorDirective set(Atom a1, Atom a2) {
        atom1 = a1;
        atom2 = a2;
        atomCount = (atom1 != null && atom2 != null) ? 2 : -1;
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