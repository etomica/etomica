package etomica.action;

import etomica.api.ISimulation;

/**
 * Interface for classes that apply some elementary action (transformation) to a
 * simulation.
 *  
 */
public interface SimulationAction extends Action {

	public ISimulation getSimulation();
}