/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import etomica.data.DataPump;
import etomica.data.DataStreamAction;

import java.util.ArrayList;

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

    public void setStreamAction(DataStreamAction action) {
        streamAction = action;
    }

    private static final long serialVersionUID = 1L;
	private final ArrayList<DataPump> dataStreamPumps;
    private DataStreamAction streamAction;
}
