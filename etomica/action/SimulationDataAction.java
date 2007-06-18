package etomica.action;

import java.util.ArrayList;

import etomica.data.DataStreamAction;

/**
 * Action that performs a call to the reset() method of a set
 * of accumulators, as specified via a list of AccumulatorManager
 * instances.
 */
public class SimulationDataAction implements Action, java.io.Serializable {

    public SimulationDataAction(DataStreamAction action) {
		dataStreamPumps = new ArrayList();
        streamAction = action;
	}
    
    public ArrayList getDataStreamPumps() {
        return dataStreamPumps;
    }

	public void actionPerformed() {
        for (int i=0; i<dataStreamPumps.size(); i++) {
            streamAction.setStart(dataStreamPumps.get(i));
            streamAction.actionPerformed();
        }
    }

    private static final long serialVersionUID = 1L;
	private final ArrayList dataStreamPumps;
    private DataStreamAction streamAction;
}
