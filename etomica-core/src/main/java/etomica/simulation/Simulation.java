/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

import java.util.HashMap;
import java.util.LinkedList;

import etomica.action.ActionIntegrate;
import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IElement;
import etomica.api.IIntegrator;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.api.ISimulationEventManager;
import etomica.api.ISpecies;
import etomica.space.ISpace;
import etomica.util.Arrays;
import etomica.util.RandomMersenneTwister;
import etomica.util.RandomNumberGeneratorUnix;

/**
 * The main class that organizes the elements of a molecular simulation.
 * Holds a single Space instance that is referenced in
 * many places to obtain spatial elements such as vectors and boundaries.
 */
public class Simulation implements ISimulation  {

    /**
     * Creates a new simulation using the given space
     */
    public Simulation(ISpace space) {
        this.space = space;
        boxList = new IBox[0];
        controller = new Controller();
        seeds = RandomNumberGeneratorUnix.getRandSeedArray();
        random = new RandomMersenneTwister(seeds);
        eventManager = new SimulationEventManager(this);
        speciesList = new ISpecies[0];
        elementSymbolHash = new HashMap<String,IElement>();
        elementAtomTypeHash = new HashMap<IElement,LinkedList<IAtomType>>();
    }

    /**
     * Returns the seeds that were used for the random number generator at
     * construction.  If the random number generator has been set manually
     * since then, this method returns null.
     */
    public int[] getRandomSeeds() {
        return seeds;
    }

    public final void addBox(IBox newBox) {
        for (int i=0; i<boxList.length; i++) {
            if (boxList[i] == newBox) {
                throw new IllegalArgumentException("Box "+newBox+" is already a part of this Simulation");
            }
        }
        boxList = (IBox[])Arrays.addObject(boxList, newBox);
        newBox.setIndex(boxList.length-1);
        for(int i=0; i<speciesList.length; i++) {
            newBox.addSpeciesNotify(speciesList[i]);
        }
        eventManager.boxAdded(newBox);
    }
    
    public final void removeBox(IBox oldBox) {
        boolean found = false;
        for (int i=0; i<boxList.length; i++) {
            if (boxList[i] == oldBox) {
                found = true;
                break;
            }
        }
        if (!found) {
            throw new IllegalArgumentException("Box "+oldBox+" is not part of this Simulation");
        }

        boxList = (IBox[])Arrays.removeObject(boxList, oldBox);

        for (int i = oldBox.getIndex(); i<boxList.length; i++) {
            boxList[i].setIndex(i);
        }
        
        eventManager.boxRemoved(oldBox);

        // notify oldBox that we no longer have it.  this will reset its index
        // to 0, so we need to do this after firing notification
        oldBox.setIndex(0);
    }
    
    /**
     * Returns an array of Boxs contained in the Simulation
     */
    public final IBox getBox(int index) {
        return boxList[index];
    }

    public int getBoxCount() {
    	return boxList.length;
    }

    /**
     * Returns the Controller used to run the simulation's Actions and 
     * Activities.
     */
	public Controller getController() {
		return controller;
	}

    /**
     * @return the space
     */
    public final ISpace getSpace() {
        return space;
    }

    public IRandom getRandom() {
        return random;
    }
    
    /**
     * Set the simulation's random number generator to the given one.
     */
    public void setRandom(IRandom newRandom) {
        seeds = null;
        random = newRandom;
    }

    public ISimulationEventManager getEventManager() {
        return eventManager;
    }

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

        for (int i=0; i<boxList.length; i++) {
            boxList[i].addSpeciesNotify(species);
        }

        // this just fires an event for listeners to receive
        eventManager.speciesAdded(species);
    }

    public void removeSpecies(ISpecies removedSpecies) {

        int index = removedSpecies.getIndex();
        
        if (speciesList[index] != removedSpecies) {
            throw new IllegalArgumentException("Species to remove not found at expected location.");
        }
        
        speciesList = (ISpecies[])Arrays.removeObject(speciesList,removedSpecies);

        for(int i = index; i < speciesList.length; i++) {
            int oldIndex = speciesList[i].getIndex();
            speciesList[i].setIndex(i);
            eventManager.speciesIndexChanged(speciesList[i], oldIndex);
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
                    eventManager.atomTypeIndexChanged(speciesList[i].getAtomType(j), oldIndex);
                }
                atomTypeMaxIndex++;
            }
        }

        for (int j = 0; j < boxList.length; j++) {
            boxList[j].removeSpeciesNotify(removedSpecies);
        }

        eventManager.speciesRemoved(removedSpecies);
        eventManager.atomTypeMaxIndexChanged(atomTypeMaxIndex);
        eventManager.speciesMaxIndexChanged(speciesList.length);
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

    /**
     * Returns the Simulation's primary integrator.  If the controller holds
     * multiple integrators, the first is returned.  If the first integrator
     * holds sub-integrators, the top-level integrator is still returned.  This
     * method assumes the controller holds an ActivityIntegrate or an
     * ActionIntegrate.
     */
    public IIntegrator getIntegrator() {
        IIntegrator integrator = null;
        IAction[] controllerActions = controller.getAllActions();
        for (int i=0; i<controllerActions.length; i++) {
            if (controllerActions[i] instanceof ActivityIntegrate) {
                integrator = ((ActivityIntegrate)controllerActions[i]).getIntegrator();
                break;
            }
            if (controllerActions[i] instanceof ActionIntegrate) {
                integrator = ((ActionIntegrate)controllerActions[i]).getIntegrator();
                break;
            }
        }
        return integrator;
    }

    private static final long serialVersionUID = 4L;
    protected final ISpace space;
    protected int[] seeds;
    protected final SimulationEventManager eventManager;
    private IBox[] boxList;
    protected IRandom random;
    private Controller controller;
    private ISpecies[] speciesList;
    private final HashMap<String,IElement> elementSymbolHash;
    private final HashMap<IElement,LinkedList<IAtomType>> elementAtomTypeHash;
}
