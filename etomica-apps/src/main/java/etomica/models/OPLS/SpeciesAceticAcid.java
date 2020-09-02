/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.OPLS;

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
 * Species for OPLS acetic acid
 * 
 * @author Hye Min Kim
 * Nov, 2011
 */
public class SpeciesAceticAcid {

    public final static int indexCH3 = 0;
    public final static int indexC = 1;
    public final static int indexDBO = 2;
    public final static int indexSBO = 3;
    public final static int indexH = 4;

    public static SpeciesGeneral create() {
        Space space = Space3D.getInstance();
        AtomType cH3Type = new AtomType(new ElementSimple("CH3", 15.0107));//mass doesn't affect anything in MC simulation
        AtomType cType = new AtomType(Carbon.INSTANCE);
        AtomType dBOType = new AtomType(Oxygen.INSTANCE, "DBO");
        AtomType sBOType = new AtomType(Oxygen.INSTANCE, "SBO");
        AtomType hType = new AtomType(Hydrogen.INSTANCE);
        return new SpeciesBuilder(space)
                .withConformation(new ConformationAceticAcid(space))
                .addCount(cH3Type, 1)
                .addCount(cType, 1)
                .addCount(dBOType, 1)
                .addCount(sBOType, 1)
                .addCount(hType, 1)
                .build();
    }
}
