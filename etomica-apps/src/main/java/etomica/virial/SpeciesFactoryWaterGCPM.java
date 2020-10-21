/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.models.water.ConformationWaterGCPM;
import etomica.models.water.SpeciesWater4P;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;


/**
 * SpeciesFactory that makes SpeciesWater
 */
public class SpeciesFactoryWaterGCPM implements SpeciesFactory, java.io.Serializable {
    public ISpecies makeSpecies(Space _space) {
        return SpeciesWater4P.create(new ConformationWaterGCPM(Space3D.getInstance()));
    }
}
