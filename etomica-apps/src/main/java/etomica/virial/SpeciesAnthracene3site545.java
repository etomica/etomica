/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.AtomType;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.Hydrogen;
import etomica.space.Space;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;

/**
 *  
 * Species Anthracene molecule
 * rigid , LJ potential, no charge, 3 site model.each benzene ring is a site
 * reference: Iwai, monte carlo sim of Naphthalene, phenathlene,anthracene in SCF 1998
 *  * @author shu
 * March, 7, 2011
 *
 */
public class SpeciesAnthracene3site545 {
    public final static int indexC1 = 0;
    public final static int indexCH1 = 1;
    public final static int indexCH2 = 2;

    public static SpeciesGeneral create(Space space) {
        //CMass is the middle site, including 4 carbons and 0 hydrogens
        //CHMass is the side site, including 5 carbons and 5 hydrogens
        double CMass = 4 * Carbon.INSTANCE.getMass();
        double CHMass = 5 * Carbon.INSTANCE.getMass() + Hydrogen.INSTANCE.getMass() * 5;
        AtomType chType = AtomType.simple("C4", CHMass);
        AtomType cType = AtomType.simple("C5H5", CMass);
        return new SpeciesBuilder(space)
                .addCount(cType, 1)
                .addCount(chType, 2)
                .withConformation(new ConformationAnthracene3site(space))
                .build();
    }
}
