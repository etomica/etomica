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
    
    public static Potential2Hard NULL = new NULL();
    public static class NULL implements Potential2Hard {
        private NULL() {}
        public double energy(AtomPair pair) {return 0.0;}
        public void bump(AtomPair pair) {}
        public double collisionTime(AtomPair pair) {return Double.MAX_VALUE;}
        public double lastCollisionVirial() {return 0.0;}
        public Space.Tensor lastCollisionVirialTensor() {return null;} //need to know D to return zero tensor
    }

    
}//end of Potential2Hard

    