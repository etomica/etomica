package etomica;

/**
 * Superclass for all Potential classes
 */
public class Hamiltonian implements java.io.Serializable {
    
    public final PotentialMaster potential;
    private final Simulation parentSimulation;
    
    public Hamiltonian(Simulation sim) {
        parentSimulation = sim;
        potential = new PotentialMaster(sim);
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