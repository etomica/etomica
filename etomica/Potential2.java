package etomica; 

/**
 * Potential acting on or within an atom, or between a pair of atoms or atom
 * groups. Contains other Potential instances that describe the specific
 * interactions between the atoms of the group(s).
 *
 * @author David Kofke
 */

 /* History of changes
  * 07/13/02 (DAK) Restructured instantiation of LRC potential
  * 07/15/02 (DAK) Constructor makes P0LRC only if instance of Potential2SoftSpherical
  * 12/06/02 (DAK) Added setIterators1A method
  * 01/27/03 (DAK) Numerous changes with redesign of Potential.
  * 07/17/03 (DAK) Made calculate method not final (so etomica.virial.P2Cluster
  * could override it)
  */

public abstract class Potential2 extends Potential {
  
    public static String VERSION = "Potential2:03.01.27/"+Potential.VERSION;
    
    protected AtomPairIterator iterator;
    
    public Potential2(SimulationElement parent) {
        this(parent, Default.TRUNCATE_POTENTIALS ? 
                        new PotentialTruncationSimple()
                      : PotentialTruncation.NULL);
      /*                  
        super(parent);
        iterator1 = new ApiIntergroup1A(parentSimulation());
        iteratorA = new ApiIntergroupAA(parentSimulation());
        if(Default.TRUNCATE_POTENTIALS) {//can't use other constructor because of "this" in constructor of PotentialTruncationSimple
            potentialTruncation = new PotentialTruncationSimple(parentSimulation().space, Default.POTENTIAL_CUTOFF_FACTOR * Default.ATOM_SIZE);
            Potential0GroupLrc lrcMaster = parentSimulation().hamiltonian.potential.lrcMaster();
            potentialTruncation.makeLrcPotential(lrcMaster, this); //adds this to lrcMaster
        } else {
            potentialTruncation = PotentialTruncation.NULL;
        }*/
    }
    public Potential2(SimulationElement parent, PotentialTruncation potentialTruncation) {
        super(2, parent, potentialTruncation);
        iterator = new Api1A(new ApiIntergroup1A(simulation()),
        						new ApiIntergroupAA(simulation()));
        if( (potentialTruncation != PotentialTruncation.NULL) && (potentialTruncation != null
            && (this instanceof Potential2SoftSpherical)) ) {
            PotentialMaster potentialMaster = simulation().hamiltonian.potential;
            potentialTruncation.makeLrcPotential(potentialMaster, this); //constructor of lrcPotential adds it to lrcMaster of potentialMaster
        }
    }
    
    public abstract double energy(AtomPair pair);
    
	/**
	 * Interface for all hard pair potentials.
	 *
	 * @author David Kofke
	 */
	public interface Hard extends Potential.Hard {
    
    	public double energy(AtomPair pair);
		/**
		 * Implements the collision dynamics.
		 * The given atoms are assumed to be at the point of collision.  This method is called
		 * to change their momentum according to the action of the collision.  Extensions can be defined to
		 * instead implement other, perhaps unphysical changes.
		 */
		public void bump(AtomPair pair);
    
		/**
		 * Computes the time of collision of the given atoms , assuming no intervening collisions.
		 * Usually assumes free-flight between collisions
		 */ 
		public double collisionTime(AtomPair pair);
    
	}//end of Potential2.Hard

	/**
	 * Methods for properties obtained for a soft, differentiable pair potential.
	 *
	 * @author David Kofke
	 */

	public interface Soft{
    
    	public double energy(AtomPair pair);
		/**
		 * Returns r dot grad(u), with any truncation applied.  Does not include
		 * division by D, to avoid repeated multiplication of this term when summing
		 * over all pairs.  Negation and division by D in most cases is required 
		 * at some point when using this quantity.
		 */
		public double virial(AtomPair pair);
    
		public double hyperVirial(AtomPair pair);
    
		public Space.Vector gradient(AtomPair pair);
    
		/**
		 * Integral used to evaluate correction to truncation of potential.
		 */
		public abstract double integral(double rC);
    
	}//end of Potential2.Soft 
	           
}//end of Potential2



