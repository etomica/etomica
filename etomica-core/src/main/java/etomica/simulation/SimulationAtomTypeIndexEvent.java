/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

import etomica.atom.AtomType;

public class SimulationAtomTypeIndexEvent extends SimulationAtomTypeEvent {

    private int index;

    public SimulationAtomTypeIndexEvent(Simulation sim, AtomType at, int index) {
        super(sim, at);
        this.index = index;
    }

    public int getIndex() {
        return index;
    }
}
