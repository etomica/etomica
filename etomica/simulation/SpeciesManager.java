package etomica.simulation;

import java.util.HashMap;
import java.util.LinkedList;

import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IEvent;
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
        elementSymbolHash = new HashMap<String,Element>();
        elementAtomTypeHash = new HashMap<Element,LinkedList<IAtomType>>();
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
            atomTypeMaxIndex += speciesList[i].getChildTypeCount();
        }
        int index = speciesList.length;
        species.setIndex(index);
        speciesList = (ISpecies[])Arrays.addObject(speciesList,species);
	    
        for(int i = 0; i < species.getChildTypeCount(); i++) {
            species.getChildType(i).setIndex(atomTypeMaxIndex++);
            atomTypeAddedNotify(species.getChildType(i));
        }

        int boxCount = sim.getBoxCount();
        for (int i=0; i<boxCount; i++) {
            sim.getBox(i).addSpeciesNotify(species);
        }

        // this just fires an event for listeners to receive
        IEvent evt = new SimulationSpeciesAddedEvent(species);
	    sim.getEventManager().fireEvent(evt);

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
            IEvent evt = new SimulationSpeciesIndexChangedEvent(speciesList[i], oldIndex);
            sim.getEventManager().fireEvent(evt);
        }

        for(int j = 0; j < removedSpecies.getChildTypeCount(); j++) {
            atomTypeRemovedNotify(removedSpecies.getChildType(j));
        }


        int atomTypeMaxIndex = 0;
        for(int i = 0; i < speciesList.length; i++) {
            for(int j = 0; j < speciesList[j].getChildTypeCount(); j++) {
                if(speciesList[i].getChildType(j).getIndex() != atomTypeMaxIndex) {
                    int oldIndex = speciesList[i].getChildType(j).getIndex();
                    speciesList[i].getChildType(j).setIndex(atomTypeMaxIndex);
                    IEvent evt = new SimulationAtomTypeIndexChangedEvent(speciesList[i].getChildType(j), oldIndex);
                    sim.getEventManager().fireEvent(evt);
                }
                atomTypeMaxIndex++;
            }
        }

        int boxCount = sim.getBoxCount();
        for (int j = 0; j < boxCount; j++) {
            sim.getBox(j).removeSpeciesNotify(removedSpecies);
        }
        IEvent evt = new SimulationSpeciesRemovedEvent(removedSpecies);
        sim.getEventManager().fireEvent(evt);

        IEvent maxEvt = new SimulationAtomTypeMaxIndexEvent(atomTypeMaxIndex);
        sim.getEventManager().fireEvent(maxEvt);

        IEvent speciesMaxEvt = new SimulationSpeciesMaxIndexEvent(speciesList.length);
        sim.getEventManager().fireEvent(speciesMaxEvt);

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
        Element newElement = newChildType.getElement();
        Element oldElement = elementSymbolHash.get(newElement.getSymbol());
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
        Element oldElement = removedType.getElement();
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
    private final HashMap<String,Element> elementSymbolHash;
    private final HashMap<Element,LinkedList<IAtomType>> elementAtomTypeHash;
    private final ISimulation sim;
}
