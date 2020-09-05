/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.atom.AtomType;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Nitrogen;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;

/**
 * 
 * 
 * Species nitrogen molecule
 * 	with 4-points charges
 * 
 * @author Tai Boon Tan
 *
 */
public class SpeciesN2 {

    public final static int indexN1 = 0;
    public final static int indexN2 = 1;
    public final static int indexP1left = 2;
    public final static int indexP2left = 3;
    public final static int indexP1right = 4;
    public final static int indexP2right = 5;

    public static SpeciesGeneral create(boolean isDynamic) {
        AtomType nType = new AtomType(Nitrogen.INSTANCE);
        AtomType pType = new AtomType(new ElementSimple("P", 1.0));
        Space space = Space3D.getInstance();
        return new SpeciesBuilder(space)
                .setDynamic(isDynamic)
                .withConformation(new ConformationNitrogen(space))
                .addCount(nType, 2)
                .addCount(pType, 4)
                .build();
    }
}
