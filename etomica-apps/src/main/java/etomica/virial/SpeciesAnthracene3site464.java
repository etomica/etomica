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
 * rigid , LJ potential, no charge, this model is very similar to 3site 545 model, only the weight of the sites are different
 * reference: Iwai, monte carlo sim of Naphthalene, phenathlene,anthracene in SCF 1998
 *  * @author shu
 * March, 7, 2011
 *
 */
public class SpeciesAnthracene3site464 {

     public static SpeciesGeneral create(Space space) {
         //CMass is the middle site, including 6 carbons and 2 hydrogens
         //CHMass is the side site, including 4 carbons and 4 hydrogens
         double CMass = 6 * Carbon.INSTANCE.getMass() + Hydrogen.INSTANCE.getMass() * 2;
         double CHMass = 4 * Carbon.INSTANCE.getMass() + Hydrogen.INSTANCE.getMass() * 4;
         AtomType chType = AtomType.simple("C4H4", CHMass);
         AtomType cType = AtomType.simple("C6H2", CMass);
         return new SpeciesBuilder(space)
                 .addCount(cType, 1)
                 .addCount(chType, 2)
                 .withConformation(new ConformationAnthracene3site(space))
                 .build();
     }
}
