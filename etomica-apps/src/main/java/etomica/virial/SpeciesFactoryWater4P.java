/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.models.water.SpeciesWater4P;
import etomica.api.ISpecies;
import etomica.config.IConformation;
import etomica.space.ISpace;


/**
 * SpeciesFactory that makes SpeciesWater
 */
public class SpeciesFactoryWater4P implements SpeciesFactory, java.io.Serializable {
    public SpeciesFactoryWater4P(IConformation conformation) {
        this.conformation = conformation;
    }
    public ISpecies makeSpecies(ISpace space) {
        SpeciesWater4P species = new SpeciesWater4P(space);
        species.setConformation(conformation);
        return species;
    }
    
    protected IConformation conformation;
}
