/*
 * History
 * Created on Aug 31, 2004 by kofke
 */
package etomica;

/**
 * Interface for hard potentials, having impulsive forces.
 */    
public interface PotentialHard {

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
            
}