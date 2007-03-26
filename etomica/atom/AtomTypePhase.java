package etomica.atom;

import etomica.simulation.SpeciesManager;


/**
 * AtomType class for the root of the AtomType tree (type of the SpeciesRoot 
 * atom)
 * @author David Kofke
 */
public class AtomTypePhase extends AtomTypeGroup {

    /**
     * Used only to create root type
     */
    public AtomTypePhase(AtomAddressManager indexManager, SpeciesManager speciesManager) {
        super(indexManager);
        this.speciesManager = speciesManager;
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
