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
        reservoirCount = -1;
    }

    /**
     * Removes the given AtomTypes associated with the given Species from the 
     * Simulation and does cleanup, including renumbering indices and firing 
     * AtomType-related event notifications.
     */
    void removeSpecies(Species removedSpecies) {
        AtomType[] agentTypes = ((AtomTypeGroup)childTypes[0]).childTypes;
        
        for (int i=0; i<agentTypes.length; i++) {
            if (agentTypes[i].getSpecies() == removedSpecies) {
                AtomTypeGroup removedType = (AtomTypeGroup)agentTypes[i];
                ((AtomTypeGroup)childTypes[0]).removeChildType(removedType);
                break;
            }
        }
        
        eventManager.fireEvent(new SimulationAtomTypeMaxIndexEvent(numChildTypes));
    }
    
    void childTypeRemovedNotify(AtomType removedType) {
        reservoirCount = 0;
        if (removedType instanceof AtomTypeLeaf) {
            Element oldElement = ((AtomTypeLeaf)removedType).getElement();
            elementSymbolHash.remove(oldElement.getSymbol());
            indexReservoir = new int[]{removedType.getIndex()};
        }
        else if (removedType instanceof AtomTypeGroup) {
            removeElements((AtomTypeGroup)removedType);
            int[] childIndices = null;
            childIndices = getChildIndices((AtomTypeGroup)removedType);
            indexReservoir = new int[childIndices.length+1];
            numChildTypes -= indexReservoir.length;
            indexReservoir[0] = removedType.getIndex();
            System.arraycopy(childIndices, 0, indexReservoir, 1, childIndices.length);
            java.util.Arrays.sort(indexReservoir);
        }
        // we can start at the agent types (our child's childTypes) because 
        // types for the SpeciesRoot and SpeciesMaster never go away
        reservoirCount = 0;
        recycleIndices(((AtomTypeGroup)childTypes[0]).childTypes);
        indexReservoir = null;
        reservoirCount = -1;
    }

    /**
     * Reassigns indices from the reservoir to the given AtomTypes.
     */
    private void recycleIndices(AtomType[] atomTypes) {
        // now iterate over remaining AtomTypes and re-use old indices that are
        // less than remaining indices
        for (int i=0; i<atomTypes.length; i++) {
            if (atomTypes[i].getIndex() >= numChildTypes) {
                int oldIndex = atomTypes[i].getIndex();
                // this triggers a call back to our requestIndex method, which will
                // return an index from the reservoir
                atomTypes[i].resetIndex();
                eventManager.fireEvent(new SimulationAtomTypeIndexChangedEvent(atomTypes[i], oldIndex));
                if (reservoirCount == indexReservoir.length) {
                    // we ran out of indices to recycle
                    return;
                }
            }
            if (atomTypes[i] instanceof AtomTypeGroup) {
                recycleIndices(((AtomTypeGroup)atomTypes[i]).childTypes);
                if (reservoirCount == indexReservoir.length) {
                    // we ran out of indices to recycle
                    return;
                }
            }
        }
    }

    /**
     * Returns an array of indices for the give parent AtomType.  The array
     * of indices does not include the given parent's index.
     */
    private static int[] getChildIndices(AtomTypeGroup atomType) {
        int[] childIndices = new int[0];
        for (int i=0; i<atomType.childTypes.length; i++) {
            AtomType childType = atomType.childTypes[i];
            if (childType instanceof AtomTypeGroup) {
                // recursion!
                int[] iChildIndices = getChildIndices((AtomTypeGroup)childType);
                int oldSize = childIndices.length;
                int newSize = oldSize + iChildIndices.length + 1;
                childIndices = Arrays.resizeArray(childIndices, newSize);
                childIndices[oldSize] = childType.getIndex();
                System.arraycopy(iChildIndices, 0, childIndices, oldSize+1, iChildIndices.length);
            }
            else {
                childIndices = Arrays.resizeArray(childIndices, childIndices.length+1);
                childIndices[childIndices.length-1] = childType.getIndex();
            }
        }
        return childIndices;
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
        if (indexReservoir == null) {
            return ++numChildTypes;
        }
        return indexReservoir[reservoirCount++];
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
    private int[] indexReservoir;
    private int reservoirCount;
    private static final long serialVersionUID = 2L;
}
