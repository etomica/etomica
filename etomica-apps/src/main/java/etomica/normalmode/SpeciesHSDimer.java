/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.AtomType;
import etomica.space.Space;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;

/**
 * 
 * 
 * Species hard-sphere dimer molecule
 * 
 * @author Tai Boon Tan
 *
 */
public class SpeciesHSDimer {

    public final static int indexAtom1 = 0;
    public final static int indexAtom2 = 1;
     public static SpeciesGeneral create(Space space, boolean isDynamic, double L) {
        return new SpeciesBuilder(space)
                .addCount(AtomType.simple("P", 1.0), 2)
                .setDynamic(isDynamic)
                .withConformation(new ConformationHSDimer(space, L))
                .build();
     }
}
