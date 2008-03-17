package etomica.api;

import etomica.action.activity.Controller;
import etomica.simulation.SimulationEventManager;
import etomica.simulation.SpeciesManager;
import etomica.space.Space;

public interface ISimulation {

    /**
     * Adds a Box to the simulation.  This method should not be called if
     * newBox is already held by the simulation.
     */
    public void addBox(IBox newBox);

    /**
     * Removes a Box to the simulation.  This method should not be called if
     * oldBox is not held by the simulation.
     */
    public void removeBox(IBox oldBox);

    /**
     * Returns an array of Boxs contained in the Simulation
     */
    public IBox[] getBoxs();

    /**
     * Returns the Controller used to run the simulation's Actions and 
     * Activities.
     */
    public IController getController();

    /**
     * Returns the Simulation's random number generator.
     */
    public IRandom getRandom();

    /**
     * Returns the Simulation's event manager, which fires events for
     * Boxs and Species being added and removed.
     */
    public SimulationEventManager getEventManager();

    /**
     * Returns the SpeciesManager, which tracks the Species in the Simulation.
     */
    public SpeciesManager getSpeciesManager();

    /**
     * @return Returns a flag indicating whether the simulation involves molecular dynamics.
     */
    public boolean isDynamic();
}