package etomica; 

/**
 * Potential acting on a single atom or atom group.
 *
 * @author David Kofke
 */
public abstract class Potential1 extends Potential {
  
    public static String VERSION = "Potential1:01.07.26/"+Potential.VERSION;
    
    protected AtomIterator iterator;
    
    public Potential1(PotentialGroup parent) {
        super(1, parent);
        iterator = simulation().iteratorFactory.makeGroupIteratorSequential();
    }
    
	public void calculate(AtomSet basis, IteratorDirective id, PotentialCalculation pc) {
		if(!enabled) return;
		switch(id.atomCount()) {
			case 0: iterator.all(basis, id, (AtomActive)pc.set(this)); break;
			case 1: singletIterator.all(basis, id, (AtomActive)pc.set(this)); break;
		}
	}//end of calculate
    
	private final AtomIteratorSinglet singletIterator = new AtomIteratorSinglet();
    
    /**
     * Convenience method for setting species.  Makes array with argument as
     * passes it to setSpecies(Species[]) method.
     * @param s The species to which this potential applies.
     */
    public void setSpecies(Species s) {
        setSpecies(new Species[] {s});
    }
    
	/**
	 * Checks that argument is allowed, and simply passes it to superclass
	 * method. No change in iterator is performed.
	 * @param species An array of exactly one element, the species to which this
	 * potential applies.
	 */
    public void setSpecies(Species[] species) {
    	if(species.length != 1) throw new IllegalArgumentException("setSpecies in Potential1 must take an array containing only one species instance");
        super.setSpecies(species);
    }

    /**
     * Returns the energy of the given atom.
     */
    public abstract double energy(Atom atom);
                      
	public void setIterator(AtomSetIterator iterator) {
		if(iterator instanceof AtomIterator) this.setIterator((AtomIterator)iterator);
		else throw new IllegalArgumentException("Inappropriate type of iterator set for potential");
	}
	public void setIterator(AtomIterator iterator) {
		this.iterator = iterator;
	}
	public AtomSetIterator getIterator() {return iterator;}
    
    /**
     * Marker interface indicating that a one-body potential is an intramolecular
     * potential, and not, e.g., a potential of interaction with an external field.
     * This is useful when computing energy changes for molecule translations and
     * rotations, for which intramolecular contributions can be ignored.
     */
    public interface Intramolecular {}
    
	/**
	 * Methods needed to describe the behavior of a hard one-body potential.  
	 * A hard potential describes impulsive interactions, in which the energy undergoes a step
	 * change at some point in the space.
	 */
 
	 /* History of changes
	  * 01/26/03 (DAK) changed to interface and put inside Potential1 class
	  * 08/26/02 (DAK) Modified calculate to handle cases of atomCount 0 or 1 in
	  * iterator directive
	  */
	public interface Hard extends Potential.Hard {

	 public double energy(Atom atom);
	 /**
	  * Implements the collision dynamics.
	  * The given atom is assumed to be at the point of collision.  This method is called
	  * to change its momentum according to the action of the collision.  Extensions can be defined to
	  * instead implement other, perhaps unphysical changes.
	  */
		public void bump(Atom atom);

	 /**
	  * Computes the time of collision of the given atom with the hard potential, assuming no intervening collisions.
	  * Usually assumes free-flight between collisions
	  */ 
		public double collisionTime(Atom atom);
            
	}  //end of Potential1.Hard

	public interface Soft {
    
    	public double energy(Atom atom);
		public Space.Vector gradient(Atom atom);
    
	}
}//end of Potential1



