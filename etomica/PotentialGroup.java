package etomica;

public interface PotentialGroup {
    
    public void addPotential(Potential p);
    
    public Simulation parentSimulation();
    
}