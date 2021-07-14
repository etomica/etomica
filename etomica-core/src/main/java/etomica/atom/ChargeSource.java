/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.atom;

import etomica.species.SpeciesManager;

public class ChargeSource {
    protected final double[] charges;

    public ChargeSource(SpeciesManager sm) {
        charges = new double[sm.getAtomTypeCount()];
    }

    public void setCharge(AtomType atomType, double q) {
        charges[atomType.getIndex()] = q;
    }

    public double getCharge(AtomType atomType) {
        return charges[atomType.getIndex()];
    }
}
