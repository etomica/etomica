/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.AlkaneEH;

import etomica.atom.AtomType;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.Hydrogen;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;

/**
 * Species for methane with explicit hydrogen
 * Bond angle = 109.5 -- not exactly, conformation is changed!!!
 * @author shu
 * 01-27-2013
 */
public class SpeciesMethane {

    public final static int indexC = 0;
    public final static int indexH1 = 1;
    public final static int indexH2 = 2;
    public final static int indexH3 = 3;
    public final static int indexH4 = 4;

     public static SpeciesGeneral create(boolean isDynamic) {
         AtomType cType = new AtomType(Carbon.INSTANCE);
         AtomType hType = new AtomType(Hydrogen.INSTANCE);
         Space space = Space3D.getInstance();
         return new SpeciesBuilder(space)
                 .setDynamic(isDynamic)
                 .withConformation(new ConformationMethane(space))
                 .addCount(cType, 1)
                 .addCount(hType, 4)
                 .build();
     }
}
