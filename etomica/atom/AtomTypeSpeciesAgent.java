package etomica.atom;

import etomica.simulation.SpeciesManager;
import etomica.species.Species;


/**
 * AtomType class for SpeciesAgents.  This lives at the top level of the
 * AtomType hierarchy and passes on notifications and requests to the
 * SpeciesManager
 * 
 * @author Andrew Schultz
 */
public class AtomTypeSpeciesAgent extends AtomTypeGroup {

    public AtomTypeSpeciesAgent(AtomAddressManager indexManager, Species species, 
            SpeciesManager speciesManager) {
        super(indexManager);
        setSpecies(species);
        this.speciesManager = speciesManager;
        index = requestIndex();
    }

    public int requestIndex() {
        return speciesManager.requestTypeIndex();
    }
    
    public void childTypeAddedNotify(AtomType childType) {
        speciesManager.atomTypeAddedNotify(childType);
    }
    
    public void childTypeRemovedNotify(AtomType childType) {
        speciesManager.atomTypeRemovedNotify(childType);
    }
    
    private static final long serialVersionUID = 1L;
    private final SpeciesManager speciesManager;
}
