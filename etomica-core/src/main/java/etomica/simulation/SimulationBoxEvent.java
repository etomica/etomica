/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

import etomica.box.Box;
import etomica.api.ISimulationBoxEvent;

public class SimulationBoxEvent extends SimulationEvent implements ISimulationBoxEvent {

    private final Box box;
    private static final long serialVersionUID = 1L;
    
    public SimulationBoxEvent(Simulation sim, Box box) {
        super(sim);
        this.box = box;
    }

    public Box getBox() {
        return box;
    }

}
