package etomica; 

/**
 * Potential acting on a single atom or atom group.
 *
 * @author David Kofke
 */
public abstract class Potential1 extends Potential {
      
    public Potential1(SimulationElement parent) {
        super(1, parent);
    }
    
    /**
     * Returns the energy of the given atom.
     */
    public abstract double energy(Atom atom);
                          
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



