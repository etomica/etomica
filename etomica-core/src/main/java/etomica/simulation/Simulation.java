/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

import etomica.action.controller.Controller;
import etomica.box.Box;
import etomica.integrator.Integrator;
import etomica.meta.annotations.IgnoreProperty;
import etomica.meta.javadoc.KeepSimJavadoc;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.species.SpeciesManager;
import etomica.util.random.IRandom;
import etomica.util.random.RandomMersenneTwister;
import etomica.util.random.RandomNumberGeneratorUnix;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * The main class that organizes the elements of a molecular simulation.
 * It contains boxes, species, an integrator, a random number generator,
 * and a controller.
 */
@KeepSimJavadoc
public class Simulation {

    protected final Space space;
    protected final SimulationEventManager eventManager;
    protected int[] seeds;
    protected IRandom random;
    private final Controller controller;

    private final List<Box> boxes;
    private SpeciesManager speciesManager;
    private final SpeciesManager.Builder smBuilder;

    /**
     * Creates a new simulation using the given space
     *
     * @param space the space used to construct Vectors etc.
     */
    public Simulation(Space space) {
        this(space, null);
    }

    public Simulation(Space space, SpeciesManager speciesManager) {
        this.space = space;
        boxes = new ArrayList<>();
        seeds = RandomNumberGeneratorUnix.getRandSeedArray();
        random = new RandomMersenneTwister(seeds);
        eventManager = new SimulationEventManager(this);
        controller = new Controller();
        this.speciesManager = speciesManager;
        this.smBuilder = speciesManager == null ? SpeciesManager.builder() : null;
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

        for(ISpecies s : getSpeciesManager().getSpeciesArray()) {
            newBox.addSpeciesNotify(s);
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
        if (speciesManager != null) {
            throw new IllegalStateException("Cannot add species after SpeciesManager is built");
        }
        this.smBuilder.addSpecies(species);
    }

    public final SpeciesManager getSpeciesManager() {
        if (speciesManager == null) {
            speciesManager = smBuilder.build();
        }
        return speciesManager;
    }

    /**
     * @return the number of Species in the Simulation.
     */
    @IgnoreProperty
    public final int getSpeciesCount() {
        return getSpeciesManager().getSpeciesCount();
    }

    public final int getAtomTypeCount() {
        return getSpeciesManager().getAtomTypeCount();
    }

    /**
     * Returns one of the simulation ISpecies. ISpecies are numbered in
     * the order that they are added, beginning with 0.
     *
     * @param index specifies the ISpecies to be returned.
     * @return the specified ISpecies.
     */
    public final ISpecies getSpecies(int index) {
        return getSpeciesManager().getSpecies(index);
    }

    public final List<ISpecies> getSpeciesList() {
        return getSpeciesManager().getSpeciesList();
    }

    /**
     * Convenience method to return the first (and often only) ISpecies
     * in the simulation.
     *
     * @return the ISpecies at index 0
     */
    public final ISpecies species() {
        return getSpeciesManager().getSpecies(0);
    }

    /**
     * Method to allow generation of unique string to identify Elements in the Simulation.
     *
     * @param symbolBase the base string.
     * @return an Element symbol starting with symbolBase that does not yet
     * exist in the Simulation.  Return values will be like "base0, base1, base2..."
     */
    public String makeUniqueElementSymbol(String symbolBase) {
        if (this.speciesManager != null) {
            throw new IllegalStateException("Too late");
        }
        return this.smBuilder.makeUniqueElementSymbol(symbolBase);
    }

    /**
     * @return the Simulation's primary integrator.  If the controller holds
     * multiple integrators, the first is returned.  If the first integrator
     * holds sub-integrators, the top-level integrator is still returned.  This
     * method assumes the controller holds an ActivityIntegrate.
     * Returns null if no integrator is found.
     */
    public Integrator getIntegrator() {
        return this.controller.getIntegrator();
    }
}
