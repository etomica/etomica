package etomica;

/**
 * Coordinator of all species agents in a phase.  
 *
 * @author David Kofke
 */
public class SpeciesMaster extends AtomGroup {
    
    private final Phase parentPhase;
    
    public SpeciesMaster(Phase p) {
        super(null, 0, AtomType.NULL); //parent, index, atomtype
//        super(null, 0, null, 0, null); //parent, index, factory, nChild, configuration
        parentPhase = p;
    }
        
    public void addSpecies(Species species) {
        SpeciesAgent agent = species.makeAgent(this, childCount+1);
        addAtom(agent);
    }
    
    public SpeciesAgent firstSpecies() {return (SpeciesAgent)firstChild();}
    public SpeciesAgent lastSpecies() {return (SpeciesAgent)lastChild();}
        
    public Phase parentPhase() {return parentPhase;}
    public Species parentSpecies() {return null;}
    public SpeciesAgent parentSpeciesAgent() {return null;}
    
    public Simulation parentSimulation() {return parentPhase.parentSimulation();}
}//end of SpeciesMaster