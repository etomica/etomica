package etomica; 

/**
 * Potential acting on a single atom or atom group.
 *
 * @author David Kofke
 */
public abstract class Potential1 extends Potential {
  
    public static String VERSION = "Potential1:01.07.26/"+Potential.VERSION;
    
    protected AtomIterator iterator;
    private Species species;
    
    public Potential1(PotentialGroup parent) {
        super(parent);
        iterator = parentSimulation().iteratorFactory.makeGroupIteratorSequential();
    }
    
	public void calculate(AtomSet basis, IteratorDirective id, PotentialCalculation pc) {
		if(!enabled) return;
		switch(id.atomCount()) {
			case 0: iterator.all(basis, id, (AtomActive)pc.set(this)); break;
			case 1: singletIterator.all(basis, id, (AtomActive)pc.set(this)); break;
		}
	}//end of calculate
    
	private final AtomIteratorSinglet singletIterator = new AtomIteratorSinglet();
    
    public void setSpecies(Species s) {
        setSpecies(new Species[] {s});
    }
    
    public void setSpecies(Species[] species) {
        switch (species.length) {
            case 1: this.species = species[0];
                    break;
            default: throw new IllegalArgumentException("Wrong number of species given in Potential1");
        }
        if(!(parentPotential() instanceof PotentialMaster)) throw new IllegalStateException("Error: Can set species only for potentials that apply at the molecule level.  Potential must have PotentialMaster as parent");
        ((PotentialMaster)parentPotential()).setSpecies(this, species);
    }
    /**
     * Returns an array of length 2 with the species to which this potential applies.
     * Returns null if no species has been set, which is the case if the potential
     * is not describing interactions between molecule-level Atoms.
     */
    public Species[] getSpecies() {
        if(species == null) return null;
        else return new Species[] {species};
    }

    /**
     * Returns the energy of the given atom.
     */
    public abstract double energy(Atom atom);
                      
    public void setIterator(AtomIterator iterator) {
        this.iterator = iterator;
    }
    public AtomIterator iterator() {return iterator;}
    
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
	public interface Hard {

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
    
		public Space.Vector gradient(Atom atom);
    
	}
}//end of Potential1



