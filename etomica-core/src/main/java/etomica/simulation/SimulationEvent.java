/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

import etomica.api.ISimulation;
import etomica.api.ISimulationEvent;


public class SimulationEvent implements ISimulationEvent, java.io.Serializable {
    
    private ISimulation simulation;
    private static final long serialVersionUID = 1L;
    
    public SimulationEvent(ISimulation sim) {
    	simulation = sim;
    }
    
    public ISimulation getSimulation() {
        return simulation;
    }
}