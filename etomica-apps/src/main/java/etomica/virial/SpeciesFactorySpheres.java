/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.chem.elements.ElementSimple;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheresMono;


/**
 * SpeciesFactory that makes SpeciesSpheresMono 
 */
public class SpeciesFactorySpheres implements SpeciesFactory, java.io.Serializable {
    public ISpecies makeSpecies(Space _space) {
        return new SpeciesSpheresMono(_space, new ElementSimple("A"));
    }
}
