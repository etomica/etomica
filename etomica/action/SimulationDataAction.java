/*
 * History
 * Created on Nov 4, 2004 by kofke
 */
package etomica.action;

import etomica.data.DataStreamAction;
import etomica.simulation.DataStreamHeader;
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
        DataStreamHeader[] streams = simulation.getDataStreams();
        for (int i=0; i<streams.length; i++) {
            Object[] clients = streams[i].getClients();
            for (int j=0; j<clients.length; j++) {
                streamAction.setStart(clients[j]);
                streamAction.actionPerformed();
            }
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
