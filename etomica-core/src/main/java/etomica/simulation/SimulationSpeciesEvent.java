/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

import etomica.api.ISimulation;
import etomica.api.ISimulationSpeciesEvent;
import etomica.api.ISpecies;

public class SimulationSpeciesEvent extends SimulationEvent
                                   implements ISimulationSpeciesEvent {

    private final ISpecies species;
    private static final long serialVersionUID = 1L;
    
    public SimulationSpeciesEvent(ISimulation sim, ISpecies species) {
        super(sim);
        this.species = species;
    }

    public ISpecies getSpecies() {
        return species;
    }

}
