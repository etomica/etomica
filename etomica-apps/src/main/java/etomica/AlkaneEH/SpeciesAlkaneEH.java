/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.AlkaneEH;


import etomica.atom.Atom;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.molecule.IMolecule;
import etomica.molecule.Molecule;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;

/**
 * Species for TraPPE-Explicit hydrogen,Siepmann
 *
 * @author shu 01-30-2013
 */
public class SpeciesAlkaneEH {
    public static SpeciesGeneral create(int numCarbons) {
        AtomType c_3Type = AtomType.simple("C3", 1.0);
        AtomType c_2Type = AtomType.simple("C2", 1.0);
        AtomType hType = AtomType.simple("H", 1.0);
        int numCH2 = numCarbons - 2;
        int numH = numCarbons * 2 + 2;
        Space space = Space3D.getInstance();
        return new SpeciesBuilder(space)
                .withConformation(new ConformationAlkaneEH(space, numCarbons))
                .addCount(c_2Type, numCH2)
                .addCount(c_3Type, 1)
                .addCount(hType, numH)
                .build();
    }
}
