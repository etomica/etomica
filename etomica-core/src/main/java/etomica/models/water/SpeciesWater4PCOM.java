/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.water;

import etomica.atom.AtomType;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;

/**
 * Species for 4-point water molecule.  This species has an explicit site for
 * the center of mass.
 */
public class SpeciesWater4PCOM {

    public final static int indexH1 = 0;
    public final static int indexH2 = 1;
    public final static int indexO = 2;
    public final static int indexM = 3;
    public final static int indexC = 4;

     public static SpeciesGeneral create(boolean isDynamic) {
        Space space = Space3D.getInstance();
         AtomType hType = new AtomType(Hydrogen.INSTANCE);
         AtomType oType = new AtomType(Oxygen.INSTANCE);
         AtomType mType = new AtomType(new ElementSimple("M", 0.0));
         AtomType cType = new AtomType(new ElementSimple("COM", 0.0));
         return new SpeciesBuilder(space)
                 .setDynamic(isDynamic)
                 .withConformation(new ConformationWaterGCPMCOM(space))
                 .addCount(hType, 2)
                 .addCount(oType, 1)
                 .addCount(mType, 1)
                 .addCount(cType, 1)
                 .build();
     }
}
