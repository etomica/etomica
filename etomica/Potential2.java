package etomica; 

/**
 * Potential acting on a pair of atoms or atom groups.  The Potential2Group
 * subclass of this potential is used to describe the inter-group interactions.
 * The Potential2Group would contain other Potential2 (or higher-body) instances
 * that describe the specific interactions between the atoms of the group.
 *
 * @author David Kofke
 */

 /* History of changes
  * 07/13/02 (DAK) Restructured instantiation of LRC potential
  * 07/15/02 (DAK) Constructor makes P0LRC only if instance of Potential2SoftSpherical
  * 12/06/02 (DAK) Added setIterators1A method
  */

public abstract class Potential2 extends Potential {
  
    public static String VERSION = "Potential2:01.07.03/"+Potential.VERSION;
    
    protected AtomPairIterator iterator;
    private Species species1, species2;
    
    public Potential2(PotentialGroup parent) {
        this(parent, Default.TRUNCATE_POTENTIALS ? 
                        new PotentialTruncationSimple(parent.parentSimulation().space)
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
    public Potential2(PotentialGroup parent, PotentialTruncation potentialTruncation) {
        super(parent, potentialTruncation);
        iterator = new Api1A(new ApiIntergroup1A(parentSimulation()),
        						new ApiIntergroupAA(parentSimulation()));
        if( (potentialTruncation != PotentialTruncation.NULL) && (potentialTruncation != null
            && (this instanceof Potential2SoftSpherical)) ) {
            PotentialMaster potentialMaster = parentSimulation().hamiltonian.potential;
            potentialTruncation.makeLrcPotential(potentialMaster, this); //constructor of lrcPotential adds it to lrcMaster of potentialMaster
        }
    }
    
    public abstract double energy(AtomPair pair);

	public final void calculate(AtomSet basis, IteratorDirective id, PotentialCalculation pc) {
	   if(!enabled) return;
	   iterator.all(basis, id, (AtomPairActive)pc.set(this));
    }//end of calculate

   public void setSpecies(Species s1, Species s2) {
    	if(s1 == s2) setSpecies(new Species[] {s1});
    	else setSpecies(new Species[] {s1, s2});
    }
    public void setSpecies(Species[] species) {
    	if(species.length == 1) {
 		   	species1 = species2 = species[0];
    	} else {
    		species1 = species[0];
    		species2 = species[1];
    	}
        if(species1 == null || species2 == null) throw new NullPointerException("Cannot set null Species in Potential2");
        if(species1 == species2) {
        	iterator = new Api1A(new ApiIntragroup1A(parentSimulation()),
									new ApiIntragroupAA(parentSimulation()));
        } else {
			iterator = new Api1A(new ApiIntergroup1A(parentSimulation()),
									new ApiIntergroupAA(parentSimulation()));
        }
		if(!(parentPotential() instanceof PotentialMaster)) throw new RuntimeException("Error: Can set species only for potentials that apply at the molecule level.  Potential must have PotentialMaster as parent");
        ((PotentialMaster)parentPotential()).setSpecies(this, species);
    }

    public void setIterator(AtomPairIterator iterator) {
        this.iterator = iterator;
    }
    public void setIterators1A(AtomPairIterator iter1, AtomPairIterator iterA) {
        iterator = new Api1A(iter1, iterA);
    }
    public AtomPairIterator iterator() {return iterator;}
    
    /**
     * Returns an array of length 2 with the species to which this potential applies.
     * Returns null if no species has been set, which is the case if the potential
     * is not describing interactions between molecule-level Atoms.
     */
    public Species[] getSpecies() {
        if(species1 == null) return null;
        else return new Species[] {species1, species2};
    }

	/**
	 * Interface for all hard pair potentials.
	 *
	 * @author David Kofke
	 */
	public interface Hard {
    
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



