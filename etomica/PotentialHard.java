package etomica;

public interface PotentialHard {
    
//public void bump(IntegratorHardAbstract.Agent agent);
    public void bump(AtomPair pair);
    
    public double lastCollisionVirial();
    
    public Space.Tensor lastCollisionVirialTensor();
    
}