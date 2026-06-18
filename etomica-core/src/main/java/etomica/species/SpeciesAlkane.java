/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.species;

import etomica.atom.AtomType;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.Hydrogen;
import etomica.config.ConformationChainZigZag;
import etomica.space.Vector;
import etomica.space3d.Space3D;

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
        return makeBuilder(numCarbons, ch3Type, ch2Type, bondL, bondTheta).build();
    }

    public static ConformationChainZigZag makeConformation() {
        return makeConformation(nominalBondL, nominalBondTheta);
    }

    public static ConformationChainZigZag makeConformation(double bondL, double bondTheta) {
        return new ConformationChainZigZag(
                Space3D.getInstance(),
                Vector.of(bondL, 0, 0),
                Vector.of(-bondL * Math.cos(bondTheta), bondL * Math.sin(bondTheta), 0)
        );
    }

    public static SpeciesBuilder makeBuilder(int numCarbons) {
        double cMass = Carbon.INSTANCE.getMass();
        double hMass = Hydrogen.INSTANCE.getMass();
        return makeBuilder(numCarbons, AtomType.simple("CH3", cMass+3*hMass),
                                       AtomType.simple("CH2", cMass+2*hMass));
    }

    public static SpeciesBuilder makeBuilder(int numCarbons, AtomType ch3Type, AtomType ch2Type) {
        return makeBuilder(numCarbons, ch3Type, ch2Type, nominalBondL, nominalBondTheta);
    }

    public static SpeciesBuilder makeBuilder(int numCarbons, AtomType ch3Type, AtomType ch2Type, double bondL, double bondTheta) {
        ConformationChainZigZag conf = makeConformation(bondL, bondTheta);
        return new SpeciesBuilder(Space3D.getInstance())
                .addCount(ch3Type, 1)
                .addCount(ch2Type, Math.max(numCarbons - 2, 0))
                .addCount(ch3Type, Math.min(1, numCarbons-1))
                .withConformation(conf);
    }
}
