package etomica;

/**
 * Coordinator of all species agents in a phase.  
 *
 * @author David Kofke
 */
public class SpeciesMaster extends AtomGroup {
    
    private final Phase parentPhase;
    private int moleculeCount, atomCount;
    
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
    
    protected void addAtomNotify(Atom atom) {
        if(atom.parentGroup() instanceof SpeciesAgent) {moleculeCount++;}
        atomCount += atom.atomCount();
    }

    protected void removeAtomNotify(Atom atom) {
        if(atom.parentGroup() instanceof SpeciesAgent) {moleculeCount--;}
        atomCount -= atom.atomCount();
    }
    
    public int atomCount() {return atomCount;}
    public int moleculeCount() {return moleculeCount;}
    
    public String signature() {return parentPhase.getName();}

}//end of SpeciesMaster