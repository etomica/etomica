/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.api.ISpecies;
import etomica.chem.elements.ElementSimple;
import etomica.space.ISpace;
import etomica.species.SpeciesSpheresRotating;


/**
 * SpeciesFactory that makes SpeciesSpheresMono 
 */
public class SpeciesFactoryOrientedSpheres implements SpeciesFactory, java.io.Serializable {
    public ISpecies makeSpecies(ISpace space) {
        return new SpeciesSpheresRotating(space, new ElementSimple("O"));
    }
}
