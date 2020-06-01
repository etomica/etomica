/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.traPPE;

import etomica.atom.AtomType;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;

/**
 * Species for methanol with satellite site (Rowley et al 2006).
 */
public class SpeciesMethanol {

     public static SpeciesGeneral create() {
        return create(false);
     }

     public static SpeciesGeneral create(boolean isDynamic) {

         double bondCH3O = 1.43; // Angstroms
         double bondOH = 0.945; // Angstroms  (Chen et al report 0.945 Angstroms..., the website says 0.95 Angstroms)
         double angleEq = 108.50*Math.PI/180; // equilibrium bond angle in radians (mcWiggle will change this appropriately)

         Vector cH3 = Vector.of(bondCH3O, 0, 0);
         Vector oxygen = Vector.of(0, 0, 0);
         Vector hydrogren = Vector.of(bondOH*Math.cos(angleEq), bondOH*Math.sin(angleEq), 0.0);

         AtomType cH3Type = new AtomType(new ElementSimple("CH3", 1.0)); // diameter taken to be CH3-CH3 equilibrium LJ distance
         AtomType oType = new AtomType(Oxygen.INSTANCE); // diameter taken to be O-O equilibrium LJ distance
         AtomType hType = new AtomType(Hydrogen.INSTANCE); // H-H equilibrium distance is not applicable

         return new SpeciesBuilder(Space3D.getInstance())
                 .addAtom(cH3Type, cH3, "CH3")
                 .addAtom(oType, oxygen, "O")
                 .addAtom(hType, hydrogren, "H")
                 .setDynamic(isDynamic)
                 .build();
     }
}
