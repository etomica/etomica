/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.atom.Atom;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomType;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Nitrogen;
import etomica.molecule.IMolecule;
import etomica.molecule.Molecule;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;

/**
 * 
 * 
 * Species nitrogen molecule (shell model) 
 * 
 * Reference: Fabianski R. et al, Calculations on the stability of low temperature solid nitrogen
 *             phases, JCP 112(15) 6745 (2000)
 *             
 * @author Tai Boon Tan
 *
 */
public class SpeciesN2ShellModel {

    public final static int indexN1 = 0;
    public final static int indexN2 = 1;
    public final static int indexCenter = 2;
    public final static int indexP1left = 3;
    public final static int indexP1right = 4;

     public static SpeciesGeneral create(boolean isDynamic) {
         AtomType nType = new AtomType(Nitrogen.INSTANCE);
         AtomType pType = new AtomType(new ElementSimple("P", 1.0));
         return new SpeciesBuilder(Space3D.getInstance())
                 .setDynamic(isDynamic)
                 .withConformation(new ConformationNitrogenShellModel(Space3D.getInstance()))
                 .addCount(nType, 2)
                 .addCount(pType, 3)
                 .build();
     }
}
