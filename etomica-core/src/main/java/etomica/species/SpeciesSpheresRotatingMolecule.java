/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.species;

import etomica.atom.Atom;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculeOriented;
import etomica.molecule.MoleculeOrientedDynamic;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Species in which molecules are made of a single atom.  The molecule itself
 * holds the orientation.
 *
 * @author Andrew Schultz
 */
public class SpeciesSpheresRotatingMolecule {

    public static SpeciesGeneral create(Space space, AtomType atomType, Vector moment, boolean isDynamic) {
        return new SpeciesBuilder(space)
                .addAtom(atomType, space.makeVector())
                .setMoleculeOriented(true)
                .setMomentOfInertia(moment)
                .setDynamic(isDynamic)
                .build();
    }
}
