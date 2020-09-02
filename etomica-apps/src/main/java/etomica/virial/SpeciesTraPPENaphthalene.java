/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.Atom;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomType;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.ElementSimple;
import etomica.molecule.IMolecule;
import etomica.molecule.Molecule;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;

/**
 *  
 * Species Naphthalene molecule
 * this is for TraPPE, the Naphthalene is rigid , LJ potential
 * reference: TraPPE 4, UA description of linear and branched alkanes and alkylbenzenes, Siepmann
 * 
 * @author shu
 * Oct, 20, 2010
 *
 */
public class SpeciesTraPPENaphthalene {

    public final static int indexC1 = 0;
    public final static int indexC2 = 1;
    public final static int indexCH1 = 2;
    public final static int indexCH2 = 3;
    public final static int indexCH3 = 4;
    public final static int indexCH4 = 5;
    public final static int indexCH5 = 6;
    public final static int indexCH6 = 7;
    public final static int indexCH7 = 8;
    public final static int indexCH8 = 9;

     public static SpeciesGeneral create(boolean isDynamic) {
         AtomType chType = new AtomType(new ElementSimple("CH", 13.0107));
         AtomType cType = new AtomType(Carbon.INSTANCE);
         Space space = Space3D.getInstance();
         return new SpeciesBuilder(space)
                 .setDynamic(isDynamic)
                 .withConformation(new ConformationNaphthaleneTraPPE(space))
                 .addCount(cType, 2)
                 .addCount(chType, 8)
                 .build();
     }
}
