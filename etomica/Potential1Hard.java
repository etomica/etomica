package etomica;
    
/**
 * Methods needed to describe the behavior of a hard one-body potential.  
 * A hard potential describes impulsive interactions, in which the energy undergoes a step
 * change at some point in the space.
 */
public abstract class Potential1Hard extends Potential1 implements PotentialHard {

    public Potential1Hard(Simulation sim) {
        super(sim);
    }
   public void bump(AtomPair pair) {bump(pair.atom1());}
   
 /**
  * Implements the collision dynamics.
  * The given atom is assumed to be at the point of collision.  This method is called
  * to change its momentum according to the action of the collision.  Extensions can be defined to
  * instead implement other, perhaps unphysical changes.
  */
    public abstract void bump(Atom atom);

 /**
  * Computes the time of collision of the given atom with the hard potential, assuming no intervening collisions.
  * Usually assumes free-flight between collisions
  */ 
    public abstract double collisionTime(Atom atom);
            
}  //end of Potential1Hard

