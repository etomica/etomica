package etomica;

/**
 * Superclass for all PotentialField and Potential classes
 */
public abstract class PotentialAbstract {
    
    public abstract double energy();
    public abstract double energy(Atom a);
    
    
//    interface Hard {
    /**
    * Implements the collision dynamics.
    * The given pair of atoms are assumed to be at the point of collision.  This method is called
    * to change their momenta according to the action of the collision.  Extensions can be defined to
    * instead implement other, perhaps unphysical changes.
    */
//        public void bump(Atom[] atoms);
    /**
    * Computes the time of collision of the given pair, assuming no intervening collisions.
    * Usually assumes free-flight between collisions
    */ 
//        public double collisionTime(Atom[] atoms);
    /**
    * Returns the collision virial from the last collision processed by this potential
    * This quantity can be used to measure the pressure
    */
//        public double lastCollisionVirial();
            
    /**
    *Returns the virial tensor from the last collision processed.  This is used to measure 
    *the pressure tensor, and eventually the surface tension
    */
//        public etomica.Space.Tensor lastCollisionVirialTensor();
//    }  //end of Potential.Hard
    
}