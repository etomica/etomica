package etomica;

public interface PotentialGroup {

    /**
     * Adds given potential to the group.  Method should not be called
     * directly, but is called when the setParentPotential method of the
     * given potential is invoked with this as the new parent.
     */
    void addPotential(Potential p);
    
    /**
     * Removes given potential from group.  Method should not be called
     * directly, but is called when the setParentPotential method of the given
     * potential is invoked with another (or null) new parent.
     */
    void removePotential(Potential p);
    
    /**
     * Returns true if this group potential contains the given potential.
     */
    public boolean contains(Potential p);
    
    public Simulation parentSimulation();
        
}