/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.species;

import etomica.atom.AtomOriented;
import etomica.atom.AtomOrientedDynamic;
import etomica.atom.AtomTypeOriented;
import etomica.atom.IAtom;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.IElement;
import etomica.simulation.Simulation;
import etomica.space.Space;

/**
 * Species in which molecules are made of a single atom of type OrientedSphere
 *
 * @author David Kofke
 * @see AtomTypeOriented
 */
public class SpeciesSpheresRotating extends SpeciesSpheresMono {

    protected boolean isAxisSymmetric = true;
    
    public SpeciesSpheresRotating(Simulation sim, Space _space) {
        this(_space, new ElementSimple(sim));
    }
    
    public SpeciesSpheresRotating(Space _space, IElement element) {
        super(_space, new AtomTypeOriented(element, _space));
    }

    public boolean isAxisSymmetric() {
        return isAxisSymmetric;
    }

    public void setAxisSymmetric(boolean isAxisSymmetric) {
        this.isAxisSymmetric = isAxisSymmetric;
    }

    protected IAtom makeLeafAtom() {
        return isDynamic ? new AtomOrientedDynamic(space, leafAtomType, isAxisSymmetric)
                         : new AtomOriented(space, leafAtomType, isAxisSymmetric);
    }
}
