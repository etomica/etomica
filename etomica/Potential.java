package etomica;

/**
 * Superclass for all Potential classes.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 8/14/02 (DAK) made parentPotential mutable, so that potential can be added/removed from a potential group
  *               added setParentPotential for this purpose.
  */
public abstract class Potential extends SimulationElement {
    
    public static String VERSION = "Potential:01.07.22";
    
    private PotentialGroup parentPotential;
	public final PotentialTruncation potentialTruncation;
	protected boolean enabled = true;
//	protected AtomSetIterator iterator;

	/**
	 * Constructor for use only by PotentialMaster subclass.
	 * @param sim Simulation instance in which potential is used.
	 */
	Potential(Simulation sim) {
		super(sim.space, Potential.class, -1);//index = -1 is arbitrary choice
		potentialTruncation = PotentialTruncation.NULL;
		if(!(this instanceof PotentialMaster)) throw new RuntimeException("Invalid attempt to instantiate potential");
	}
	
    public Potential(PotentialGroup parent) {
    	this(parent, PotentialTruncation.NULL);
    }
    public Potential(PotentialGroup parent, PotentialTruncation truncation) {
        super(parent.parentSimulation(), Potential.class);
        potentialTruncation = truncation;
        parentPotential = parent;
        parentPotential.addPotential(this);
    }

    public final PotentialGroup parentPotential() {return parentPotential;}
     
    /**
     * Adds this potential to the given potential group.  No action is taken
     * if the new parent is the same as the current parent.  Otherwise, if
     * current parent is not null, this is removed from it.  Then this is
     * added via addPotential to the new parent.  If new parent is null, this
     * is removed from current parent and no new parent is set.
     */
    public void setParentPotential(PotentialGroup newParent) {
        if(newParent == parentPotential) return;
        if(parentPotential != null) parentPotential.removePotential(this);
        parentPotential = newParent;
        if(newParent != null) parentPotential.addPotential(this);
    }
    
    public abstract void calculate(AtomSet basis, IteratorDirective id, PotentialCalculation pc);
            
//	public void calculate(AtomSet basis, IteratorDirective id, PotentialCalculation pc) {
//		if(!enabled) return;
//		iterator.all(basis, id, pc);
//	}
	
    public void setSpecies(Species[] s) {
    	if(parentPotential instanceof PotentialMaster) {
    		((PotentialMaster)parentPotential).setSpecies(this, s);
    	} else {
    		throw new RuntimeException("Error: attempt to set species for potential that is not at top level");
    	}
    }
    
	/**
	 * Returns the enabled flag, which if false will cause potential to not
	 * contribute to any potential calculations.
	 * @return boolean The current value of the enabled flag
	 */
	public boolean isEnabled() {
		return enabled;
	}

	/**
	 * Sets the enabled flag, which if false will cause potential to not
	 * contribute to any potential calculations.
	 * @param enabled The enabled value to set
	 */
	public void setEnabled(boolean enabled) {
		this.enabled = enabled;
	}
	
	/**
	 * Accessor method for potential cutoff implementation.
	 */
	public PotentialTruncation getTruncation() {return potentialTruncation;}

	/**
	 * Marker interface for Null potentials, which are defined to have no action.
	 */
	public interface Null {}
        

}//end of Potential