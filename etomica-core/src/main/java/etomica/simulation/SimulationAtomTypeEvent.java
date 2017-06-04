/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

import etomica.atom.IAtomType;
import etomica.api.ISimulationAtomTypeEvent;

public class SimulationAtomTypeEvent extends SimulationEvent
                            implements ISimulationAtomTypeEvent {
    
    private final IAtomType atomType;
    private static final long serialVersionUID = 1L;
    
    public SimulationAtomTypeEvent(Simulation sim, IAtomType atomType) {
        super(sim);
        this.atomType = atomType;
    }

    public IAtomType getAtomType() {
        return atomType;
    }

}
