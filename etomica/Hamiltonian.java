package etomica;

/**
 * Superclass for all Potential classes
 */
public class Hamiltonian implements java.io.Serializable {
    
    public final PotentialMaster potential;
    private final Simulation parentSimulation;
    
    //may need to have a different potential hierarchy for each integrator,
    //if multiple integrators are to be run on separate threads
    
    public Hamiltonian(Simulation sim, PotentialMaster potential) {
        parentSimulation = sim;
        this.potential = potential;
    }

    /**
     * Resets integrators in all phases in which the potential has an agent.
     */
/*    public void resetIntegrators() {
        Iterator e = agents.values().iterator();
        while(e.hasNext()) {
            PotentialAgent agent = (PotentialAgent)e.next();
            if(agent.parentPhase().integrator() != null)
                agent.parentPhase().integrator().reset();
        }
    }*/
    
}//end of Potential