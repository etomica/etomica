package etomica;

/**
 * Superclass for all Potential classes, which define how the atoms in the
 * system interact with each other.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 06/16/03 (DAK) Revised to permit SimulationElement in constructor.
  * 01/27/03 (DAK) Large set of changes in revision of design of Potential
  * 08/14/02 (DAK) made parentPotential mutable, so that potential can be
  * added/removed from a potential group; added setParentPotential for this
  * purpose.
  */
public abstract class Potential extends SimulationElement {
    
    private PotentialGroup parentPotential;
	public final PotentialTruncation potentialTruncation;
	protected boolean enabled = true;
	public final int nBody;
	private Potential0Lrc p0Lrc;

	/**
	 * Constructor for use only by PotentialMaster subclass.
	 * @param sim Simulation instance in which potential is used.
	 */
	Potential(Simulation sim) {
		super(sim, Potential.class);
		nBody = 0;
		potentialTruncation = PotentialTruncation.NULL;
		if(!(this instanceof PotentialMaster)) throw new RuntimeException("Invalid attempt to instantiate potential");
	}
		
	/**
	 * Constructor with default potential truncation given
	 * as PotentialTruncation. NULL.
	 * @param nBody number of atoms to which potential is applied at a time
	 * @param parent simulation element (usually a PotentialGroup) in which this
	 * potential resides
	 */
    public Potential(int nBody, SimulationElement parent) {
    	this(nBody, parent, PotentialTruncation.NULL);
    }
    /**
     * General constructor for a potential instance
     * @param nBody number of atoms to which this potential applies at a time;
     * for example with a pair potential nBody = 2; for a single-body potential,
     * nBody = 1.
     * @param parent potential group in which this potential reside
     * @param truncation instance of a truncation class that specifies the
     * scheme for truncating the potential
     */
    public Potential(int nBody, SimulationElement parent, PotentialTruncation truncation) {
        super(parent, Potential.class);
        this.nBody = nBody;
        potentialTruncation = truncation;
        if(parent instanceof PotentialGroup) {
	        parentPotential = (PotentialGroup)parent;
	        parentPotential.addPotential(this);
        }
    }

    public abstract double energy(Atom[] atoms);
    
    public abstract void setPhase(Phase phase);
    
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
 	/**
	 * Accessor method for potential cutoff implementation.
	 */
	public PotentialTruncation getTruncation() {return potentialTruncation;}

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
	}//end of interface Hard
	
	/**
	 * Returns the zero-body potential used to apply a long-range correction
	 * for truncation of this potential.
	 * @return Potential0Lrc
	 */
	public Potential0Lrc getLrc() {
		return p0Lrc;
	}

	/**
	 * Sets the zero-body potential used to apply a long-range correction (lrc)
	 * for truncation of this potential.  This is invoked in the constructor of
	 * Potential0Lrc, which itself is routinely invoked during the construction
	 * of this potential, so this method is declared final to guard against
	 * subclasses performing some action that is inappropriate while this class
	 * is being constructed.
	 * @param p0Lrc The lrc potential to set
	 */
	final void setLrc(Potential0Lrc p0Lrc) {
		this.p0Lrc = p0Lrc;
	}

}//end of Potential

