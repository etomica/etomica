package etomica;

public interface PotentialHard {
    
//public void bump(IntegratorHardAbstract.Agent agent);
    public void bump(AtomPair pair);
    
    public double lastCollisionVirial();
    
    public Space.Tensor lastCollisionVirialTensor();
    
    /**
     * Instance of hard pair potential corresponding to no interaction between atoms.
     */
    public static PotentialHard NULL = new NULL();
    static class NULL implements PotentialHard, Potential.Null {
        private NULL() {}
        public void bump(AtomPair pair) {}
        public double lastCollisionVirial() {return 0.0;}
        public Space.Tensor lastCollisionVirialTensor() {return null;} //need to know D to return zero tensor
    }
}