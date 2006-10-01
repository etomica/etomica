package etomica.atom;

import etomica.simulation.SimulationAtomTypeAddedEvent;
import etomica.simulation.SimulationEventManager;
import etomica.species.Species;
import etomica.util.Arrays;


/**
 * AtomType class for the root of the AtomType tree (type of the SpeciesRoot 
 * atom)
 * @author David Kofke
 */
public class AtomTypeRoot extends AtomTypeGroup {

    /**
     * Used only to create root type
     */
    AtomTypeRoot(AtomAddressManager indexManager) {
        super(indexManager);
        numChildTypes = 0;
    }

    void setSpeciesRoot(SpeciesRoot newRoot) {
        if (root != null) {
            throw new IllegalStateException("SpeciesRoot can only be called from the SpeciesRoot constructor");
        }
        root = newRoot;
    }
    
    /**
     * Removes the given AtomTypes associated with the given Species from the 
     * Simulation and does cleanup, including renumbering indices and firing 
     * AtomType-related event notifications.
     */
    void removeSpecies(Species removedSpecies) {
        AtomType[] agentTypes = ((AtomTypeGroup)childTypes[0]).childTypes;
        AtomType removedType = null;
        int startIndex = -1;  // index of first AtomType removed
        int stopIndex = -1;   // index of first AtomType not removed
        for (int i=0; i<agentTypes.length; i++) {
            if (startIndex != -1) {
                // we've already removed the Species.  renumber the rest of the indices
                // the indices should be requested in the same order they had been
                agentTypes[i].resetIndex();
            }
            if (agentTypes[i].getSpecies() == removedSpecies) {
                removedType = agentTypes[i];
                startIndex = removedType.getIndex();
                ((AtomTypeGroup)childTypes[0]).childTypes = (AtomType[])Arrays.removeObject(((AtomTypeGroup)childTypes[0]).childTypes,removedType);
                agentTypes = ((AtomTypeGroup)childTypes[0]).childTypes;
                if (i < agentTypes.length) {
                    // if it was the last newTypeIndex is -1, which is OK
                    stopIndex = ((AtomTypeGroup)childTypes[0]).childTypes[i].getIndex();
                }
                // reset numChildTypes to the number that come before the removed type.
                // we'll renumber the ones after it.
                numChildTypes = startIndex - 1;
                i--;
            }
        }
        
        eventManager.fireEvent(new SimulationAtomTypeCompactedEvent(startIndex, stopIndex));
    }
    
    int requestIndex() {
        return ++numChildTypes;
    }
    
    public void setEventManager(SimulationEventManager newEventManager) {
        if (eventManager != null) {
            throw new RuntimeException("EventManager should only be set once");
        }
        eventManager = newEventManager;
    }
    
    protected void childTypeAddedNotify(AtomType newChildType) {
        if (eventManager != null) {
            // could be null if this is the SpeciesRoot type
            eventManager.fireEvent(new SimulationAtomTypeAddedEvent(newChildType));
        }
    }
    
    private SpeciesRoot root;
    private SimulationEventManager eventManager;
    private int numChildTypes;
    private static final long serialVersionUID = 1L;
}
