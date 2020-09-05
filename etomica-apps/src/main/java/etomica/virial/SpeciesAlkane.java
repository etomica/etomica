/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.AtomType;
import etomica.config.ConformationChainZigZag2;
import etomica.molecule.IMolecule;
import etomica.molecule.Molecule;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.species.SpeciesSpheresHetero;

public class SpeciesAlkane {

    // parameters for TraPPE
    public static final double nominalBondL = 1.54;
    public static final double nominalBondTheta = Math.PI * 114 / 180;

    public static SpeciesGeneral create(int numCarbons) {
        return create(numCarbons, AtomType.simple("CH3", 15), AtomType.simple("CH2", 14));
    }

    public static SpeciesGeneral create(int numCarbons, AtomType ch3Type, AtomType ch2Type) {
        return create(numCarbons, ch3Type, ch2Type, nominalBondL, nominalBondTheta);
    }

    public static SpeciesGeneral create(int numCarbons, AtomType ch3Type, AtomType ch2Type, double bondL, double bondTheta) {
        Space space = Space3D.getInstance();
        ConformationChainZigZag2 conf = new ConformationChainZigZag2(
                space,
                Vector.of(bondL, 0, 0),
                Vector.of(-bondL * Math.cos(bondTheta), bondL * Math.sin(bondTheta), 0)
        );
        return new SpeciesBuilder(Space3D.getInstance())
                .addCount(ch3Type, Math.min(2, numCarbons))
                .addCount(ch2Type, Math.max(numCarbons - 2, 0))
                .withConformation(conf)
                .build();
    }
}
