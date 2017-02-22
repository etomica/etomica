/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

import etomica.api.ISimulation;
import etomica.api.ISimulationSpeciesIndexEvent;
import etomica.api.ISpecies;

public class SimulationSpeciesIndexEvent extends SimulationSpeciesEvent
                                 implements ISimulationSpeciesIndexEvent {

    private int index;
    
    public SimulationSpeciesIndexEvent(ISimulation sim, ISpecies species,
                                       int index) {
        super(sim, species);
        this.index = index;
    }
    
    public int getIndex() {
        return index;
    }
}
