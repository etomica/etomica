package etomica;

public interface Potential2Hard extends PotentialHard {
    
    
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
    
}//end of Potential2Hard

    