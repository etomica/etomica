package etomica.action;

import etomica.simulation.Simulation;

/**
 * Interface for classes that apply some elementary action (transformation) to a
 * simulation.
 *  
 */
public interface SimulationAction extends Action {

	public void setSimulation(Simulation sim);

	public Simulation getSimulation();
}