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
  */

public abstract class Potential2 extends Potential {
  
    public static String VERSION = "Potential2:03.01.27/"+Potential.VERSION;
    
    protected AtomPairIterator iterator;
    
    public Potential2(PotentialGroup parent) {
        this(parent, Default.TRUNCATE_POTENTIALS ? 
                        new PotentialTruncationSimple(parent.simulation().space)
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

	public final void calculate(AtomSet basis, IteratorDirective id, PotentialCalculation pc) {
	   if(!enabled) return;
	   if(Default.EXPLICIT_LOOP) {
		switch(basis.nBody()) {
		   case 1: iterator.setBasis((Atom)basis, (Atom)basis); break;
		   case 2: iterator.setBasis(((AtomPair)basis).atom1(), ((AtomPair)basis).atom2());
		}
	   	   iterator.reset(id);
	   	   pc.set(this).doCalculation(iterator, this);
	   } else {
		   iterator.all(basis, id, (AtomPairActive)pc.set(this));
	   }
    }//end of calculate

	/**
	 * Convenience method for setting species.  Simply forms array with
	 * arguments and passes it to setSpecies(Species[]) method.
	 */
   public void setSpecies(Species s1, Species s2) {
    	if(s1 == s2) setSpecies(s1);
    	else setSpecies(new Species[] {s1, s2});
    }
    
	/**
	 * Convenience method for setting species.  Simply forms array with
	 * arguments and passes it to setSpecies(Species[]) method.
	 */
   public void setSpecies(Species s) {
    	setSpecies(new Species[] {s});
    }
    
    /**
     * Extends superclass setSpecies method to instantiate iterators based on
     * number of species in array.  If only one species is given, potential is
     * assumed to be intra-molecular, and ApiIntraGroup iterators are used; if
     * two species are given, potential is assumed to be inter-molecular, and
     * ApiInterGroup iterators are used.
     * @see etomica.Potential#setSpecies(Species[])
     */
    public void setSpecies(Species[] species) {
        if(species.length < 2) {  
        	super.setSpecies(species);      	
        	iterator = new Api1A(new ApiIntragroup1A(simulation()),
									new ApiIntragroupAA(simulation()));
        } else if(species[0] == species[1]) {
        	super.setSpecies(new Species[] {species[0]});
			iterator = new Api1A(new ApiIntragroup1A(simulation()),
									new ApiIntragroupAA(simulation()));        		
        } else {
			super.setSpecies(species);
			iterator = new Api1A(new ApiIntergroup1A(simulation()),
									new ApiIntergroupAA(simulation()));
        }
    }

	public void setIterator(AtomSetIterator iterator) {
		if(iterator instanceof AtomPairIterator) this.setIterator((AtomPairIterator)iterator);
		else throw new IllegalArgumentException("Inappropriate type of iterator set for potential");
	}
    public void setIterator(AtomPairIterator iterator) {
        this.iterator = iterator;
    }
    public AtomSetIterator getIterator() {return iterator;}
    
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
	
	//static instance of this class is made in Potential.Hard
	//inconsistency is indicated if this class is defined in Potential instead of here
	final static class HardNull implements Potential1.Hard, Potential2.Hard, Potential.Null {
		public void bump(Atom a) {}
		public double collisionTime(Atom a) {return Double.MAX_VALUE;}
		public double energy(Atom a) {return 0.0;}
		public void bump(AtomPair pair) {}
		public double collisionTime(AtomPair pair) {return Double.MAX_VALUE;}
		public double energy(AtomPair pair) {return 0.0;}
		public double lastCollisionVirial() {return 0.0;}
		public Space.Tensor lastCollisionVirialTensor() {throw new etomica.exception.MethodNotImplementedException();} //need to know D to return zero tensor
	}//end of NULL
           
}//end of Potential2



