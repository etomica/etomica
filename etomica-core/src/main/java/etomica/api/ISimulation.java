/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.api;

import etomica.box.Box;
import etomica.integrator.Integrator;
import etomica.util.random.IRandom;

/**
 * The ISimulation is an interface for a simulation, containing boxes, species,
 * a random number generator and an integrator.
 * 
 * @author Andrew Schultz
 */
public interface ISimulation {

    /**
     * Adds a Box to the simulation.  This method should not be called if
     * newBox is already held by the simulation.
     */
    public void addBox(Box newBox);

    /**
     * Removes a Box to the simulation.  This method should not be called if
     * oldBox is not held by the simulation.
     */
    public void removeBox(Box oldBox);

    /**
     * Returns number of boxes contained in the Simulation
     */
    public int getBoxCount();

    /**
     * Returns Box specified by index contained in the Simulation
     */
    public Box getBox(int index);

    /**
     * Adds species to the list of all ISpecies in the simulation, and
     * adds notifies all IBoxes of the new ISpecies.
     */
    public void addSpecies(ISpecies species);

    /**
     * Removes the given ISpecies from the ISimulation.
     */
    public void removeSpecies(ISpecies removedSpecies);

    /**
     * Returns the number of Species in the Simulation.
     */
    public int getSpeciesCount();

    /**
     * Returns the Species in the Simulation for the specified index.
     */
    public ISpecies getSpecies(int index);

    /**
     * Returns the integrator for this simulation.
     */
    public Integrator getIntegrator();

    /**
     * Returns the Simulation's random number generator.
     */
    public IRandom getRandom();

    /**
     * Returns the Simulation's event manager, which fires events for
     * Boxes and Species being added and removed.
     */
    public ISimulationEventManager getEventManager();
}
