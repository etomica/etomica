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
    
	public final PotentialTruncation potentialTruncation;
	private final int nBody;
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
    }

    public abstract double energy(Atom[] atoms);
    
    public abstract void setPhase(Phase phase);
    
    public final int nBody() {return nBody;}
    
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
		  * Implements the collision dynamics.
		  * The given atom(s) is assumed to be at the point of collision.  This method is called
		  * to change their momentum according to the action of the collision.  Extensions can be defined to
		  * instead implement other, perhaps unphysical changes.
		  */
			public void bump(Atom[] atom);

		 /**
		  * Computes the time of collision of the given atom(s) with the hard potential, assuming no intervening collisions.
		  * Usually assumes free-flight between collisions.
		  */ 
			public double collisionTime(Atom[] atom);
	            
	}//end of interface Hard

	/**
	 * Methods for properties obtained for a soft, differentiable pair potential.
	 *
	 * @author David Kofke
	 */
	public interface Soft {
		   
    	public double energy(Atom[] atoms);
		/**
		 * Returns r dot grad(u), with any truncation applied.  Does not include
		 * division by D, to avoid repeated multiplication of this term when summing
		 * over all pairs.  Negation and division by D in most cases is required 
		 * at some point when using this quantity.
		 */
		public Space.Vector gradient(Atom[] pair);
    
	}
	
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

