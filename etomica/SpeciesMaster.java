package etomica;

/**
 * Coordinator of all species agents in a phase.  
 *
 * @author David Kofke
 */
public class SpeciesMaster extends AtomGroup {
    
    private final Phase parentPhase;
    
    public SpeciesMaster(Phase p) {
        super(null, 0); //parent, index
//        super(null, 0, null, 0, null); //parent, index, factory, nChild, configuration
        parentPhase = p;
    }
        
    public void addSpecies(Species species) {
        Species.Agent agent = species.makeAgent(this, childCount++, null, 0, null);
        addAtom(agent);
    }
        
    public Phase parentPhase() {return parentPhase;}
    
    public Simulation parentSimulation() {return parentPhase.parentSimulation();}
}//end of SpeciesMaster