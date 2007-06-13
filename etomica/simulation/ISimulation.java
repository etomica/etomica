package etomica.simulation;

import etomica.action.activity.Controller;
import etomica.phase.Phase;
import etomica.util.IRandom;

public interface ISimulation {

    /**
     * Adds a Phase to the simulation.  This method should not be called if
     * newPhase is already held by the simulation.
     */
    public void addPhase(Phase newPhase);

    /**
     * Removes a Phase to the simulation.  This method should not be called if
     * oldPhase is not held by the simulation.
     */
    public void removePhase(Phase oldPhase);

    /**
     * Returns an array of Phases contained in the Simulation
     */
    public Phase[] getPhases();

    /**
     * Returns the Controller used to run the simulation's Actions and 
     * Activities.
     */
    public Controller getController();

    /**
     * @return Returns a flag indicating whether the simulation involves
     * molecular dynamics.
     */
    public boolean isDynamic();

    /**
     * Returns the Simatulion's random number generator.
     */
    public IRandom getRandom();

    /**
     * Returns the Simulation's event manager, which fires events for
     * Phases and Species being added and removed.
     */
    public SimulationEventManager getEventManager();

    /**
     * Returns the SpeciesManager, which tracks the Species in the Simulation.
     */
    public SpeciesManager getSpeciesManager();
}