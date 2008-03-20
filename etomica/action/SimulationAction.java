package etomica.action;

import etomica.api.IAction;
import etomica.api.ISimulation;

/**
 * Interface for classes that apply some elementary action (transformation) to a
 * simulation.
 *  
 */
public interface SimulationAction extends IAction {

	public ISimulation getSimulation();
}