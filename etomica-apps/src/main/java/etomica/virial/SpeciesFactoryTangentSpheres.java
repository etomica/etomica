/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.chem.elements.ElementSimple;
import etomica.config.IConformation;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheres;

/**
 * SpeciesFactory that makes a tangent sphere species.
 */
public class SpeciesFactoryTangentSpheres implements SpeciesFactory, java.io.Serializable {
    public SpeciesFactoryTangentSpheres(int nA, IConformation conformation) {
        this.nA = nA;
        this.conformation = conformation;
    }
    
    public ISpecies makeSpecies(Space space) {
        return new SpeciesSpheres(nA, new ElementSimple("TS"), conformation, space);
    }
    
    private static final long serialVersionUID = 1L;
    private final int nA;
    private final IConformation conformation;
}
