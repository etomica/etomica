/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.rowley;


import etomica.atom.AtomType;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;

/**
 * Species for methanol with satellite site (Rowley et al 2006).
 */
public class SpeciesMethanol {

    public final static int indexO = 0;
    public final static int indexaC = 1;
    public final static int indexaH = 2; // ahType
    public final static int indexH1 = 3; // hType
    public final static int indexH2a = 4; // hType
    public final static int indexH2b = 5; // hType
    public final static int indexX = 6;
    public static SpeciesGeneral create(boolean pointCharges) {
        AtomType oType = new AtomType(Oxygen.INSTANCE);
        AtomType acType = new AtomType(Carbon.INSTANCE);
        AtomType ahType = new AtomType(Hydrogen.INSTANCE, "AH");
        AtomType hType = new AtomType(Hydrogen.INSTANCE);
        AtomType xType = new AtomType(new ElementSimple("X", 1.0));
        Space space = Space3D.getInstance();
        return new SpeciesBuilder(space)
                .withConformation(new ConformationMethanol(space, pointCharges))
                .addCount(oType, 1)
                .addCount(acType, 1)
                .addCount(ahType, 1)
                .addCount(hType, 4)
                .addCount(xType, 1)
                .build();
    }

}
