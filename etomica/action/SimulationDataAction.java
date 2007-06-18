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

    public SimulationDataAction(Simulation sim, DataStreamAction action) {
		simulation = sim;
        streamAction = action;
	}

	public void actionPerformed() {
        DataStreamHeader[] streams = simulation.getDataStreams();
        for (int i=0; i<streams.length; i++) {
            streamAction.setStart(streams[i].getClient());
            streamAction.actionPerformed();
        }
    }

    private static final long serialVersionUID = 1L;
	private final Simulation simulation;
    private DataStreamAction streamAction;
}
