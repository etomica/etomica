package etomica.action.controller;

/**
 * A process conducted as part of a simulation. Typically consists of repeated integration steps.
 * An Activity is given to {@link Controller} as part of a Simulation setup.
 */
public abstract class Activity {

    public void restart() { }

    public abstract void runActivity(Controller.ControllerHandle handle);

}
