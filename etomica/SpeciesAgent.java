package etomica;

/**
 * The SpeciesAgent is a representative of the species in each phase.
 * The agent handles addition, deletion, link-list ordering, counting, etc. of 
 * molecules in a phase.  Each phase has an agent from every species instance.
 * 
 * @author David Kofke
 */

public class SpeciesAgent extends AtomGroup {

    private final Species parentSpecies;
    
    public SpeciesAgent(Species s, SpeciesMaster parent, int index, int nMolecules,
                        Configuration initializer) {
        super(parent, index, AtomType.NULL, s.moleculeFactory(), nMolecules, initializer, true);
        parentSpecies = s;
    }
        
    public final Species parentSpecies() {return parentSpecies;}
    
    public final SpeciesAgent parentSpeciesAgent() {return this;}
    
    public SpeciesAgent nextSpecies() {return (SpeciesAgent)nextAtom();}
    public int moleculeCount() {return childCount();}
    public void setNMolecules(int n) {setNAtoms(n);}
    
    public final int depth() {return 1;}
        
} //end of SpeciesAgent
