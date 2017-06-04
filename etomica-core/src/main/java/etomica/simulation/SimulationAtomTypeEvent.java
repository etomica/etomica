/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

import etomica.api.ISimulationAtomTypeEvent;
import etomica.atom.AtomType;

public class SimulationAtomTypeEvent extends SimulationEvent
                            implements ISimulationAtomTypeEvent {
    
    private static final long serialVersionUID = 1L;
    private final AtomType atomType;

    public SimulationAtomTypeEvent(Simulation sim, AtomType atomType) {
        super(sim);
        this.atomType = atomType;
    }

    public AtomType getAtomType() {
        return atomType;
    }

}
