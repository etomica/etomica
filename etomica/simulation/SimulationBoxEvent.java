/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

import etomica.api.IBox;
import etomica.api.ISimulation;
import etomica.api.ISimulationBoxEvent;

public class SimulationBoxEvent extends SimulationEvent implements ISimulationBoxEvent {

    private final IBox box;
    private static final long serialVersionUID = 1L;
    
    public SimulationBoxEvent(ISimulation sim, IBox box) {
        super(sim);
        this.box = box;
    }

    public IBox getBox() {
        return box;
    }

}
