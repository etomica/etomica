/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.water;

import etomica.atom.Atom;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomType;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.molecule.IMolecule;
import etomica.molecule.Molecule;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;

/**
 * Species for 3-point water molecule.
 */
public class SpeciesWater3P {
    public static SpeciesGeneral create() {
        return create(false);
    }

    public static SpeciesGeneral create(boolean isDynamic) {
        AtomType hType = new AtomType(Hydrogen.INSTANCE);
        AtomType oType = new AtomType(Oxygen.INSTANCE);
        return new SpeciesBuilder(Space3D.getInstance())
                .withConformation(new ConformationWater3P(Space3D.getInstance()))
                .addCount(hType, 2)
                .addCount(oType, 1)
                .setDynamic(isDynamic)
                .build();
    }
}
