/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.AtomType;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.Hydrogen;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;

/**
 *
 * Species Phenanthrene molecule
 * rigid , LJ potential, no charge, 3 site model.each benzene ring is a site
 * reference: Iwai, monte carlo sim of Naphthalene, phenathlene,anthracene in SCF 1998
 * @author shu
 * March, 9, 2011
 *
 */
public class SpeciesPh3site545 {

    public final static int indexC1 = 0;
    public final static int indexCH1 = 1;
    public final static int indexCH2 = 2;

    public static SpeciesGeneral create(boolean isDynamic) {
        Space space = Space3D.getInstance();
        double CMass = 4 * Carbon.INSTANCE.getMass();
        double CHMass = 5 * Carbon.INSTANCE.getMass() + Hydrogen.INSTANCE.getMass() * 5;
        AtomType chType = AtomType.simple("C5H5", CHMass);
        AtomType cType = AtomType.simple("C4", CMass);
        return new SpeciesBuilder(space)
                .setDynamic(isDynamic)
                .withConformation(new ConformationPh3site(space))
                .addCount(cType, 1)
                .addCount(chType, 2)
                .build();

    }
}
