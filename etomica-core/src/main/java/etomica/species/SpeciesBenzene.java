/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.species;

import etomica.atom.AtomType;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.Hydrogen;
import etomica.config.ConformationBenzene;
import etomica.space3d.Space3D;

public class SpeciesBenzene {

    // parameters for TraPPE
    public static final double nominalBondL = 1.4;

    public static SpeciesGeneral create() {
        return create(AtomType.simple("CH", Carbon.INSTANCE.getMass()+Hydrogen.INSTANCE.getMass()));
    }

    public static SpeciesGeneral create(AtomType chType) {
        return create(chType, nominalBondL);
    }

    public static SpeciesGeneral create(AtomType chType, double bondL) {
        return makeBuilder(chType, bondL).build();
    }

    public static ConformationBenzene makeConformation() {
        return makeConformation(nominalBondL);
    }

    public static ConformationBenzene makeConformation(double bondL) {
        return new ConformationBenzene(bondL);
    }

    public static SpeciesBuilder makeBuilder() {
        double cMass = Carbon.INSTANCE.getMass();
        double hMass = Hydrogen.INSTANCE.getMass();
        return makeBuilder(AtomType.simple("CH", cMass+hMass));
    }

    public static SpeciesBuilder makeBuilder(AtomType chType) {
        return makeBuilder(chType, nominalBondL);
    }

    public static SpeciesBuilder makeBuilder(AtomType chType, double bondL) {
        ConformationBenzene conf = makeConformation(bondL);
        return new SpeciesBuilder(Space3D.getInstance())
                .addCount(chType, 6)
                .withConformation(conf);
    }
}
