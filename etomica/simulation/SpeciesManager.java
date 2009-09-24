package etomica.simulation;

import java.util.HashMap;
import java.util.LinkedList;

import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IElement;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.api.ISpeciesManager;
import etomica.chem.elements.Element;
import etomica.util.Arrays;

/**
 * The SpeciesManager manages Species and AtomTypes on behalf of the
 * Simulation.
 * 
 * @author Andrew Schultz
 */
public class SpeciesManager implements java.io.Serializable, ISpeciesManager {

    public SpeciesManager(ISimulation sim) {
        this.sim = sim;
        speciesList = new ISpecies[0];
        elementSymbolHash = new HashMap<String,IElement>();
        elementAtomTypeHash = new HashMap<IElement,LinkedList<IAtomType>>();
    }

    /* (non-Javadoc)
	 * @see etomica.simulation.ISpeciesManager#addSpecies(etomica.api.ISpecies)
	 */
    public void addSpecies(ISpecies species) {

        int atomTypeMaxIndex = 0;

        for (int i=0; i<speciesList.length; i++) {
            if (speciesList[i] == species) {
                throw new IllegalArgumentException("Species already exists");
            }
            atomTypeMaxIndex += speciesList[i].getAtomTypeCount();
        }
        int index = speciesList.length;
        species.setIndex(index);
        speciesList = (ISpecies[])Arrays.addObject(speciesList,species);
	    
        for(int i = 0; i < species.getAtomTypeCount(); i++) {
            species.getAtomType(i).setIndex(atomTypeMaxIndex++);
            atomTypeAddedNotify(species.getAtomType(i));
        }

        int boxCount = sim.getBoxCount();
        for (int i=0; i<boxCount; i++) {
            sim.getBox(i).addSpeciesNotify(species);
        }

        // this just fires an event for listeners to receive
        ((SimulationEventManager)sim.getEventManager()).speciesAdded(species);

    }

    /* (non-Javadoc)
	 * @see etomica.simulation.ISpeciesManager#removeSpecies(etomica.api.ISpecies)
	 */
    public void removeSpecies(ISpecies removedSpecies) {

        int index = removedSpecies.getIndex();
        
        if (speciesList[index] != removedSpecies) {
            throw new IllegalArgumentException("Species to remove not found at expected location.");
        }
        
        speciesList = (ISpecies[])Arrays.removeObject(speciesList,removedSpecies);

        for(int i = index; i < speciesList.length; i++) {
            int oldIndex = speciesList[i].getIndex();
            speciesList[i].setIndex(i);
            ((SimulationEventManager)sim.getEventManager()).speciesIndexChanged(speciesList[i], oldIndex);
        }

        for(int j = 0; j < removedSpecies.getAtomTypeCount(); j++) {
            atomTypeRemovedNotify(removedSpecies.getAtomType(j));
        }


        int atomTypeMaxIndex = 0;
        for(int i = 0; i < speciesList.length; i++) {
            for(int j = 0; j < speciesList[j].getAtomTypeCount(); j++) {
                if(speciesList[i].getAtomType(j).getIndex() != atomTypeMaxIndex) {
                    int oldIndex = speciesList[i].getAtomType(j).getIndex();
                    speciesList[i].getAtomType(j).setIndex(atomTypeMaxIndex);
                    ((SimulationEventManager)sim.getEventManager()).atomTypeIndexChanged(speciesList[i].getAtomType(j), oldIndex);
                }
                atomTypeMaxIndex++;
            }
        }

        int boxCount = sim.getBoxCount();
        for (int j = 0; j < boxCount; j++) {
            sim.getBox(j).removeSpeciesNotify(removedSpecies);
        }
        ((SimulationEventManager)sim.getEventManager()).speciesRemoved(removedSpecies);

        ((SimulationEventManager)sim.getEventManager()).atomTypeMaxIndexChanged(atomTypeMaxIndex);

        ((SimulationEventManager)sim.getEventManager()).speciesMaxIndexChanged(speciesList.length);

    }

    /* (non-Javadoc)
	 * @see etomica.simulation.ISpeciesManager#boxAddedNotify(etomica.api.IBox)
	 */
    public void boxAddedNotify(IBox newBox) {
        for(int i=0; i<speciesList.length; i++) {
            newBox.addSpeciesNotify(speciesList[i]);
        }
    }

    public int getSpeciesCount() {
    	return speciesList.length;
    }

    public ISpecies getSpecies(int index) {
        return speciesList[index];
    }

    protected void atomTypeAddedNotify(IAtomType newChildType) {
        IElement newElement = newChildType.getElement();
        IElement oldElement = elementSymbolHash.get(newElement.getSymbol());
        if (oldElement != null && oldElement != newElement) {
            // having two AtomTypes with the same Element is OK, but having
            // two Elements with the same symbol is not allowed.
            throw new IllegalStateException("Element symbol "+newElement.getSymbol()+" already exists in this simulation as a different element");
        }
        // remember the element so we can check for future duplication
        elementSymbolHash.put(newElement.getSymbol(), newElement);
        LinkedList<IAtomType> atomTypeList = elementAtomTypeHash.get(newElement);
        if (atomTypeList == null) {
            atomTypeList = new LinkedList<IAtomType>();
            elementAtomTypeHash.put(newElement, atomTypeList);
        }
        atomTypeList.add(newChildType);
    }

    protected void atomTypeRemovedNotify(IAtomType removedType) {
        // remove the type's element from our hash 
        IElement oldElement = removedType.getElement();
        elementSymbolHash.remove(oldElement.getSymbol());
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
        // this will get replaced by the actual Element when it gets added via childTypeAddedNotify
        elementSymbolHash.put(symbolBase+n, null);
        return symbolBase+n;
    }

    private static final long serialVersionUID = 1L;
    private ISpecies[] speciesList;
//    private IAtomTypeLeaf[] atomTypeList;
    private final HashMap<String,IElement> elementSymbolHash;
    private final HashMap<IElement,LinkedList<IAtomType>> elementAtomTypeHash;
    private final ISimulation sim;
}
