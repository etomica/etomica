/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.AtomType;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.ElementSimple;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;

/**
 *
 * Species Phenanthrene molecule
 * this is for TraPPE, the Phenanthrene is rigid , LJ potential
 * reference: TraPPE 4, UA description of linear and branched alkanes and alkylbenzenes, Siepmann
 * modified from species of phenanthrene C14H10
 * @author shu
 * March.19.2011
 *
 */
public class SpeciesTraPPEPhenanthrene {

    public final static int indexC1 = 0;
    public final static int indexC2 = 1;
    public final static int indexC3 = 2;
    public final static int indexC4 = 3;
    public final static int indexCH1 = 4;
    public final static int indexCH2 = 5;
    public final static int indexCH3 = 6;
    public final static int indexCH4 = 7;
    public final static int indexCH5 = 8;
    public final static int indexCH6 = 9;
    public final static int indexCH7 = 10;
    public final static int indexCH8 = 11;
    public final static int indexCH9 = 12;
    public final static int indexCH10 = 13;

    public static SpeciesGeneral create(boolean isDynamic) {
        AtomType chType = new AtomType(new ElementSimple("CH", 13.0107));
        AtomType cType = new AtomType(Carbon.INSTANCE);
        Space space = Space3D.getInstance();
        return new SpeciesBuilder(space)
                .setDynamic(isDynamic)
                .withConformation(new ConformationPhenanthreneTraPPE(space))
                .addCount(cType, 4)
                .addCount(chType, 10)
                .build();
    }
}
