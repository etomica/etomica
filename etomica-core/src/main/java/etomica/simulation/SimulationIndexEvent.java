/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

public class SimulationIndexEvent extends SimulationEvent {

    private int index;

    public SimulationIndexEvent(Simulation sim, int index) {
        super(sim);
        this.index = index;
    }

    public int getIndex() {
        return index;
    }
}
