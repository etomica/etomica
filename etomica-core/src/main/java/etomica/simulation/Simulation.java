/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.IElement;
import etomica.integrator.Integrator;
import etomica.meta.annotations.IgnoreProperty;
import etomica.meta.javadoc.KeepSimJavadoc;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.util.logging.SimLogger;
import etomica.util.random.IRandom;
import etomica.util.random.RandomMersenneTwister;
import etomica.util.random.RandomNumberGeneratorUnix;

import java.util.*;

/**
 * The main class that organizes the elements of a molecular simulation.
 * It contains boxes, species, an integrator, a random number generator,
 * and a controller.
 */
@KeepSimJavadoc
public class Simulation {
    protected final SimLogger LOG = SimLogger.create();

    protected final Space space;
    protected final SimulationEventManager eventManager;
    private final Map<String, IElement> elementSymbolHash;
    private final Map<IElement, LinkedList<AtomType>> elementAtomTypeHash;
    protected int[] seeds;
    protected IRandom random;
    private final Controller controller;

    private final List<ISpecies> speciesList;
    private final List<Box> boxes;

    /**
     * Creates a new simulation using the given space
     *
     * @param space the space used to construct Vectors etc.
     */
    public Simulation(Space space) {
        LOG.info("Initializing Simulation");
        this.space = space;
        boxes = new ArrayList<>();
        controller = new Controller();
        seeds = RandomNumberGeneratorUnix.getRandSeedArray();
        random = new RandomMersenneTwister(seeds);
        eventManager = new SimulationEventManager(this);
        speciesList = new ArrayList<>();
        elementSymbolHash = new HashMap<>();
        elementAtomTypeHash = new HashMap<>();
    }

    /**
     * @return the seeds that were used for the random number generator at
     * construction.  If the random number generator has been set manually
     * since then, this method returns null.
     */
    public int[] getRandomSeeds() {
        return seeds;
    }

    /**
     * Get all boxes in the simulation.
     * <p>
     * Attempts to add or remove boxes from this list will throw UnsupportedOperationException, you must
     * use the addBox and removeBox methods on Simulation.
     *
     * @return a read-only list of Boxes
     */
    public final List<Box> getBoxes() {
        return Collections.unmodifiableList(this.boxes);
    }

    /**
     * Convenience method to retrieve the first (and typically only) Box in the simulation.
     *
     * @return the Box at index 0
     */
    public final Box box() {
        return this.boxes.get(0);
    }

    /**
     * Adds a Box to the simulation.  This method may not be called until all Species have
     * been added to the simulation.
     *
     * @param newBox the Box being added.
     * @throws IllegalArgumentException if newBox was already added to the simulation.
     */
    public final Box addBox(Box newBox) {
        if (boxes.contains(newBox)) {
            throw new IllegalArgumentException("Box " + newBox + " is already a part of this Simulation");
        }
        boxes.add(newBox);
        newBox.setIndex(boxes.size() - 1);

        for(ISpecies aSpeciesList : speciesList) {
            newBox.addSpeciesNotify(aSpeciesList);
        }
        eventManager.boxAdded(newBox);
        return newBox;
    }

    /**
     * Creates a new Box with the default Boundary and adds it to the Simulation.
     *
     * @return the new Box.
     */
    public Box makeBox() {
        Box box = new Box(space);
        this.addBox(box);
        return box;
    }

    /**
     * Creates a new Box and adds it to the Simulation.
     *
     * @param boundary the boundary to use when constructing the Box.
     * @return the new Box.
     */
    public Box makeBox(Boundary boundary) {
        Box box = new Box(boundary, space);
        this.addBox(box);
        return box;
    }

    /**
     * Removes a Box to the simulation.
     *
     * @param oldBox the Box being removed.
     * @throws IllegalArgumentException if oldBox was not previously added to the simulation.
     */
    public final void removeBox(Box oldBox) {
        boolean found = boxes.remove(oldBox);
        if(found) {
            for (int i = oldBox.getIndex(); i < boxes.size(); i++) {
                boxes.get(i).setIndex(i);
            }

            eventManager.boxRemoved(oldBox);

            // notify oldBox that we no longer have it.  this will reset its index
            // to 0, so we need to do this after firing notification
            oldBox.setIndex(0);
        } else {
            throw new IllegalArgumentException("Box " + oldBox + " is not part of this Simulation");
        }

    }

    /**
     * Returns one of the simulation Boxes. Boxes are numbered in
     * the order that they are added, beginning with 0.
     *
     * @param index specifies the Box to be returned.
     * @return the specified Box.
     */
    public final Box getBox(int index) {
        return boxes.get(index);
    }

    /**
     * @return the number of boxes that have been added to the Simulation.
     */
    @IgnoreProperty
    public final int getBoxCount() {
        return boxes.size();
    }

    /**
     * @return the Controller used to run the simulation's Actions and
     * Activities.
     */
    @IgnoreProperty
    public final Controller getController() {
        return controller;
    }

    /**
     * @return the Space that was specified in the constructor.
     */
    public final Space getSpace() {
        return space;
    }

    /**
     * @return the Simulation's random number generator.
     */
    public final IRandom getRandom() {
        return random;
    }

    /**
     * Set the simulation's random number generator to the given one.
     *
     * @param newRandom the new random number generator.
     */
    public void setRandom(IRandom newRandom) {
        seeds = null;
        random = newRandom;
    }

    /**
     * @return the Simulation's event manager, which fires events for
     * Boxes and Species being added and removed.
     */
    @IgnoreProperty
    public final SimulationEventManager getEventManager() {
        return eventManager;
    }

    /**
     * Adds species to the list of all ISpecies in the simulation, and
     * notifies all Boxes of the addition.
     *
     * @param species the Species being added.
     * @throws IllegalArgumentException if species was already added.
     * @throws IllegalStateException if a Box has already been added to the simulation.
     */
    public final void addSpecies(ISpecies species) {
        if (boxes.size() > 0) {
            throw new IllegalStateException("Cannot add species after adding a box");
        }
        if(speciesList.contains(species)) {
            throw new IllegalArgumentException("Species already exists");
        }

        int atomTypeMaxIndex = 0;

        for (ISpecies s : speciesList) {
            atomTypeMaxIndex += s.getAtomTypeCount();
        }

        int index = speciesList.size();
        species.setIndex(index);
        speciesList.add(species);

        for (AtomType atomType : species.getAtomTypes()) {
            atomType.setIndex(atomTypeMaxIndex++);
            atomTypeAddedNotify(atomType);
        }

        for(Box box : boxes) {
            box.addSpeciesNotify(species);
        }
    }

    /**
     * @return the number of Species in the Simulation.
     */
    @IgnoreProperty
    public final int getSpeciesCount() {
        return speciesList.size();
    }

    /**
     * Returns one of the simulation ISpecies. ISpecies are numbered in
     * the order that they are added, beginning with 0.
     *
     * @param index specifies the ISpecies to be returned.
     * @return the specified ISpecies.
     */
    public final ISpecies getSpecies(int index) {
        return speciesList.get(index);
    }

    public final List<ISpecies> getSpeciesList() {
        return Collections.unmodifiableList(speciesList);
    }

    /**
     * Convenience method to return the first (and often only) ISpecies
     * in the simulation.
     *
     * @return the ISpecies at index 0
     */
    public final ISpecies species() {
        return speciesList.get(0);
    }

    private void atomTypeAddedNotify(AtomType newChildType) {
        IElement newElement = newChildType.getElement();
        IElement oldElement = elementSymbolHash.get(newElement.getSymbol());
        if (oldElement != null && oldElement != newElement) {
            // having two AtomTypes with the same Element is OK, but having
            // two Elements with the same symbol is not allowed.
            throw new IllegalStateException("Element symbol " + newElement.getSymbol() + " already exists in this simulation as a different element");
        }
        // remember the element so we can check for future duplication
        elementSymbolHash.put(newElement.getSymbol(), newElement);
        LinkedList<AtomType> atomTypeList = elementAtomTypeHash.get(newElement);
        if (atomTypeList == null) {
            atomTypeList = new LinkedList<AtomType>();
            elementAtomTypeHash.put(newElement, atomTypeList);
        }
        atomTypeList.add(newChildType);
    }

    /**
     * Method to allow generation of unique string to identify Elements in the Simulation.
     *
     * @param symbolBase the base string.
     * @return an Element symbol starting with symbolBase that does not yet
     * exist in the Simulation.  Return values will be like "base0, base1, base2..."
     */
    public String makeUniqueElementSymbol(String symbolBase) {
        int n = 0;
        while (elementSymbolHash.containsKey(symbolBase + n)) {
            n++;
        }
        // reserve this symbol so future calls to makeUniqueElementSymbol won't return it
        // this will get replaced by the actual Element when it gets added via childTypeAddedNotify
        elementSymbolHash.put(symbolBase + n, null);
        return symbolBase + n;
    }

    /**
     * @return the Simulation's primary integrator.  If the controller holds
     * multiple integrators, the first is returned.  If the first integrator
     * holds sub-integrators, the top-level integrator is still returned.  This
     * method assumes the controller holds an ActivityIntegrate.
     * Returns null if no integrator is found.
     */
    public Integrator getIntegrator() {
        Integrator integrator = null;
        IAction[] controllerActions = controller.getAllActions();
        for (IAction controllerAction : controllerActions) {
            if (controllerAction instanceof ActivityIntegrate) {
                integrator = ((ActivityIntegrate) controllerAction).getIntegrator();
                break;
            }
        }
        return integrator;
    }
}
