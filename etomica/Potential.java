package etomica;

/**
 * Superclass for all Potential classes, which define how the atoms in the
 * system interact with each other.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 01/27/03 (DAK) Large set of changes in revision of design of Potential
  * 08/14/02 (DAK) made parentPotential mutable, so that potential can be
  * added/removed from a potential group added setParentPotential for this
  * purpose.
  */
public abstract class Potential extends SimulationElement {
    
    public static String VERSION = "Potential:03.01.27";
    
    private PotentialGroup parentPotential;
	public final PotentialTruncation potentialTruncation;
	protected boolean enabled = true;
	private Species[] species;
	public final int nBody;
//	protected AtomSetIterator iterator;

	/**
	 * Constructor for use only by PotentialMaster subclass.
	 * @param sim Simulation instance in which potential is used.
	 */
	Potential(Simulation sim) {
		super(sim.space, Potential.class, -1);//index = -1 is arbitrary choice
		nBody = 0;
		potentialTruncation = PotentialTruncation.NULL;
		if(!(this instanceof PotentialMaster)) throw new RuntimeException("Invalid attempt to instantiate potential");
	}
	
    public Potential(int nBody, PotentialGroup parent) {
    	this(nBody, parent, PotentialTruncation.NULL);
    }
    public Potential(int nBody, PotentialGroup parent, PotentialTruncation truncation) {
        super(parent.parentSimulation(), Potential.class);
        this.nBody = nBody;
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
	
	/**
	 * Sets the species to which this potential applies, if it is being used for
	 * molecule-level interactions.  Subclasses may extend this method
	 * to perform actions to instantiate appropriate iterators; but overriding
	 * methods must at some point invoke this superclass method to make sure
	 * species is registered appropriately with potential master.  An exception
	 * is thrown if the parent of this potential is not the potential master
	 * (only child potentials of potential master apply to the molecule- level
	 * interactions).
	 */
    public void setSpecies(Species[] s) {
    	if(s == null || s.length == 0) throw new IllegalArgumentException("Error: setSpecies called without specifying any species instances");
    	if(s.length > nBody) throw new IllegalArgumentException("Error:  Attempting to associate potential with more species than can be defined for it");
    	if(parentPotential instanceof PotentialMaster) {
    		((PotentialMaster)parentPotential).setSpecies(this, s);
    	} else {
    		throw new RuntimeException("Error: Can set species only for potentials that apply at the molecule level.  Potential must have PotentialMaster as parent");
    	}
    	species = new Species[s.length];
    	System.arraycopy(s, 0, species, 0, s.length);
    }
    /**
     * Returns the species to which this potential applies, if it is
     * defined for molecule-level interactions.  Returns null if no species has
     * been set, which is the case if the potential is not describing
     * interactions between molecule-level Atoms.
     */
    public Species[] getSpecies() {return species;}
    
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
	 * Sets the iterator that defines the atoms to which this potential
	 * applies.  Iterator should be appropriate type for the concrete subclass
	 * of this class.
	 */
	public abstract void setIterator(AtomSetIterator iterator);
	
	/**
	 * Accessor method for the iterator that defines the atoms to which this
	 * potential applies.
	 */
	public abstract AtomSetIterator getIterator();
	
	/**
	 * Marker interface for Null potentials, which are defined to have no action.
	 */
	public interface Null {}
    
    /**
	 * Interface for hard potentials, having impulsive forces.
     */    
	public interface Hard {
    
    	/**
    	 * Value of the virial from the most recent collision.
    	 * @return double virial value
    	 */
		public double lastCollisionVirial();
    
    	/**
    	 * Value of the virial from the most recent collision, decomposed into
    	 * it tensoral elements.
    	 * @return Tensor
    	 */
		public Space.Tensor lastCollisionVirialTensor();
    
		/**
		 * Instance of hard pair potential corresponding to no interaction between atoms.
		 */
		public static Hard NULL = new NULL();
		static class NULL implements Potential.Hard, Potential.Null {
			private NULL() {}
			public double lastCollisionVirial() {return 0.0;}
			public Space.Tensor lastCollisionVirialTensor() {throw new etomica.exception.MethodNotImplementedException();} //need to know D to return zero tensor
		}//end of NULL
	}//end of interface Hard
}//end of Potential