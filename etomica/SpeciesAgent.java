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
    
    public SpeciesAgent(Species s, SpeciesMaster parent, int index, int nChild,
                        CoordinateInitializer initializer) {
        super(parent, index, s.atomFactory(), nChild, initializer, true);
        parentSpecies = s;
    }
        
    public final Species parentSpecies() {return parentSpecies;}
    
    public final byte depth() {return 1;}
        
} //end of SpeciesAgent
