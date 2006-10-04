package etomica.atom;

import java.util.HashMap;
import java.util.LinkedList;

import etomica.chem.elements.Element;
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
        elementSymbolHash = new HashMap();
        elementAtomTypeHash = new HashMap();
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
                removeElements((AtomTypeGroup)removedType);
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

    /**
     * Removes all elements from the element hash which are children of the
     * given parent AtomType
     */
    private void removeElements(AtomTypeGroup oldParentType) {
        AtomType[] oldChildTypes = oldParentType.getChildTypes();
        for (int i=0; i<childTypes.length; i++) {
            if (oldChildTypes[i] instanceof AtomTypeLeaf) {
                Element oldElement = ((AtomTypeLeaf)oldChildTypes[i]).getElement();
                elementSymbolHash.remove(oldElement.getSymbol());
            }
            else if (oldChildTypes[i] instanceof AtomTypeGroup) {
                removeElements((AtomTypeGroup)oldChildTypes[i]);
            }
        }
    }
    
    int requestIndex() {
        return ++numChildTypes;
    }
    
    /**
     * Sets the event manager for the AtomType.  This should only be called by 
     * SpeciesRoot (which actually owns the event manager).
     */
    public void setEventManager(SimulationEventManager newEventManager) {
        if (eventManager != null) {
            throw new RuntimeException("EventManager should only be set once");
        }
        eventManager = newEventManager;
    }
    
    protected void childTypeAddedNotify(AtomType newChildType) {
        if (newChildType instanceof AtomTypeLeaf) {
            Element newElement = ((AtomTypeLeaf)newChildType).getElement();
            Element oldElement = (Element)elementSymbolHash.get(newElement.getSymbol());
            if (oldElement != null && oldElement != newElement) {
                // having two AtomTypes with the same Element is OK, but having
                // two Elements with the same symbol is not allowed.
                throw new IllegalStateException("Element symbol "+newElement.getSymbol()+" already exists in this simulation as a different element");
            }
            elementSymbolHash.put(newElement.getSymbol(), newElement);
            LinkedList atomTypeList = (LinkedList)elementAtomTypeHash.get(newElement);
            if (atomTypeList == null) {
                atomTypeList = new LinkedList();
                elementAtomTypeHash.put(newElement, atomTypeList);
            }
            atomTypeList.add(newChildType);
        }
        if (eventManager != null) {
            // could be null if this is the SpeciesRoot type
            eventManager.fireEvent(new SimulationAtomTypeAddedEvent(newChildType));
        }
    }
    
    /**
     * Returns an Element symbol starting with symbolBase that does not yet 
     * exist in the Simulation.  Return values will be like "base0, base1, base2..." 
     */
    public String makeUniqueElementSymbol(String symbolBase) {
        int n = 0;
        while (elementSymbolHash.containsKey(symbolBase+n)) {
            n++;
        }
        // reserve this symbol so future calls to makeUniqueElementSymbol won't return it
        // this will get repalced by the actual Element when it gets added via childTypeAddedNotify
        elementSymbolHash.put(symbolBase+n, null);
        return symbolBase+n;
    }
    
    private final HashMap elementSymbolHash;
    private final HashMap elementAtomTypeHash;
    private SimulationEventManager eventManager;
    private int numChildTypes;
    private static final long serialVersionUID = 2L;
}
