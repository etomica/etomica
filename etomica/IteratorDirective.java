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
    
    private Direction direction;
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
        setDirection(direction);
    }
    public IteratorDirective(Direction direction, Atom atom) {
    	this(direction);
    	targetAtoms = atom;
    }
    
    /**
     * Puts directive in default state of no atoms specified, up direction, no
     * potential criteria applied, no LRC included.
     */
    public void clear() {
        setDirection(UP);
        targetAtoms = null;
        potentialCriteriaHead = null;
        includeLrc = false;
    }
    
    public final void setDirection(Direction direction) {
        this.direction = direction;
    }

    public final Direction direction() {return direction;}
    
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
    
    public void setTargetAtoms(AtomSet atoms) {
        if (Debug.ON && atoms == null) throw new IllegalArgumentException("Target atoms array must not be null");
        targetAtoms = atoms;
    }
    public AtomSet getTargetAtoms() {return targetAtoms;}
    
    public AtomSet targetAtoms;
    
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