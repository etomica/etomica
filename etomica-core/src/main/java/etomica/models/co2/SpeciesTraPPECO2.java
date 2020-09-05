/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.co2;

import etomica.atom.AtomType;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.Oxygen;
import etomica.space.Space;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;

/**
 * 
 * 
 * Species CO2 molecule
 *this is for TraPPE, the CO2 is rigid , LJ potential + QQ 
 * 
 * @author Shu Yang
 * Oct, 20, 2010
 *
 */
public class SpeciesTraPPECO2 {
    public final static int indexC = 0;
    public final static int indexOleft = 1;
    public final static int indexOright = 2;

    public static SpeciesGeneral create(Space space) {
         return create(space, false);
     }

    public static SpeciesGeneral create(Space space, boolean isDynamic) {
         AtomType cType = new AtomType(Carbon.INSTANCE);
         AtomType oType = new AtomType(Oxygen.INSTANCE);
         return new SpeciesBuilder(space)
                 .addCount(cType, 1)
                 .addCount(oType, 2)
                 .setDynamic(isDynamic)
                 .withConformation(new ConformationCO2(space))
                 .build();
     }
}
