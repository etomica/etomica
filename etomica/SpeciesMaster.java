package etomica;

/**
 * Coordinator of all species agents in a phase.  
 *
 * @author David Kofke
 */
public class SpeciesMaster extends AtomGroup {
    
    private final Phase parentPhase;
    
    public SpeciesMaster(Phase p) {
        super(p.parentSimulation().space(), null, AtomType.NULL); //space, parent, atomtype
        parentPhase = p;
    }
        
    public void addSpecies(Species species) {
        SpeciesAgent agent = species.makeAgent(this);
        addAtom(agent);
        parentPhase.getConfiguration().initializeCoordinates(this);
    }
    
    public SpeciesAgent firstSpecies() {return (SpeciesAgent)firstChild();}
    public SpeciesAgent lastSpecies() {return (SpeciesAgent)lastChild();}
        
    public Phase parentPhase() {return parentPhase;}
    public Species parentSpecies() {return null;} //never called
    public SpeciesAgent parentSpeciesAgent() {return null;}//never called
    public Simulation parentSimulation() {return parentPhase.parentSimulation();}
}//end of SpeciesMaster