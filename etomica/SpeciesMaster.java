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
        updateCounts();
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

    /**
     * Updates all the molecule and atoms counts for this phase.
     * Does not include wall atoms in the count.
     */
    public void updateCounts() {
        moleculeCount = 0;
        atomCount = 0;
        for(SpeciesAgent s=firstSpecies(); s!=null; s=s.nextSpecies()) {
            if(s.firstLeafAtom().type instanceof AtomType.Wall) continue;
            moleculeCount++;
        }
        for(Atom a=firstLeafAtom(); a!=null; a=a.nextAtom()) {
            if(a.type instanceof AtomType.Wall) continue;
            atomCount++;
        }
    }//end of updateCounts

}//end of SpeciesMaster