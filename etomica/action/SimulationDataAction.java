/*
 * History
 * Created on Nov 4, 2004 by kofke
 */
package etomica.action;

import etomica.data.DataStreamAction;
import etomica.data.DataStuff;
import etomica.simulation.Simulation;

/**
 * Action that performs a call to the reset() method of a set
 * of accumulators, as specified via a list of AccumulatorManager
 * instances.
 */
public class SimulationDataAction implements Action, java.io.Serializable {

	/**
	 * 
	 */
	public SimulationDataAction(Simulation sim, DataStreamAction action) {
		simulation = sim;
        streamAction = action;
	}

	public void actionPerformed() {
        DataStuff[] stuffs = simulation.getDataStreams();
        for (int i=0; i<stuffs.length; i++) {
            streamAction.setStart(stuffs[i]);
        }
    }

	public String getLabel() {
		return label;
	}

	public void setLabel(String label) {
		this.label = label;
	}

	private final Simulation simulation;
    private DataStreamAction streamAction;
	private String label = "Simulation Data Action";
}
