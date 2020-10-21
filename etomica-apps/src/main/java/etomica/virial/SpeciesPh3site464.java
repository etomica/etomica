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
 * rigid , LJ potential, no charge, 3 site model.each benzene ring is a site, 464 model
 * reference: Iwai, monte carlo sim of Naphthalene, phenathlene,anthracene in SCF 1998
 * @author shu
 * March, 9, 2011
 *
 */
public class SpeciesPh3site464 {

    public final static int indexC1 = 0;
    public final static int indexCH1 = 1;
    public final static int indexCH2 = 2;

     public static SpeciesGeneral create(boolean isDynamic) {
        Space space = Space3D.getInstance();
         // the virial simulation  doesnt care this , but still I added this info here FYI
         //CMass is the middle site, including 6 carbons and 2 hydrogens
         //CHMass is the side site, including 4 carbons and 4 hydrogens
         double CMass = 6 * Carbon.INSTANCE.getMass() + Hydrogen.INSTANCE.getMass() * 2;
         double CHMass = 4 * Carbon.INSTANCE.getMass() + Hydrogen.INSTANCE.getMass() * 4;
         AtomType chType = AtomType.simple("C4H4", CHMass);
         AtomType cType = AtomType.simple("C6H2", CMass);
         return new SpeciesBuilder(space)
                 .setDynamic(isDynamic)
                 .withConformation(new ConformationPh3site(space))
                 .addCount(cType, 1)
                 .addCount(chType, 2)
                 .build();

     }
}
