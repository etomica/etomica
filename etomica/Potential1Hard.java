package etomica;
    
/**
 * Methods needed to describe the behavior of a hard one-body potential.  
 * A hard potential describes impulsive interactions, in which the energy undergoes a step
 * change at some point in the space.
 */
public interface Potential1Hard extends PotentialHard {

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

