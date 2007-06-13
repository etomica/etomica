package etomica.simulation;

import etomica.action.activity.Controller;
import etomica.phase.Phase;
import etomica.util.IRandom;

public interface ISimulation {

    public void addPhase(Phase newPhase);

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
     * @return Returns a flag indicating whether the simulation involves molecular dynamics.
     */
    public boolean isDynamic();

    public IRandom getRandom();

    public SimulationEventManager getEventManager();

}