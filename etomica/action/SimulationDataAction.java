package etomica.action;

import java.util.ArrayList;

import etomica.data.DataPump;
import etomica.data.DataStreamAction;

/**
 * Action that performs a call to the reset() method of a set
 * of accumulators, as specified via a list of AccumulatorManager
 * instances.
 */
public class SimulationDataAction implements IAction, java.io.Serializable {

    public SimulationDataAction(DataStreamAction action) {
		dataStreamPumps = new ArrayList<DataPump>();
        streamAction = action;
	}
    
    public ArrayList<DataPump> getDataStreamPumps() {
        return dataStreamPumps;
    }

	public void actionPerformed() {
        for (int i=0; i<dataStreamPumps.size(); i++) {
            streamAction.setStart(dataStreamPumps.get(i));
            streamAction.actionPerformed();
        }
    }

    private static final long serialVersionUID = 1L;
	private final ArrayList<DataPump> dataStreamPumps;
    private DataStreamAction streamAction;
}
