/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.freeenergy.npath;

import etomica.atom.IAtomList;
import etomica.potential.IPotentialAtomic;
import etomica.potential.PotentialCalculation;

/**
 * Potential calculation whose sole purpose in life is to invoke getDUDW
 *
 * Created by andrew on 5/8/17.
 */
public class PotentialCalculationDUDW implements PotentialCalculation {

    protected double sum = 0;

    public void reset() {
        sum = 0;
    }

    public double getSum() {
        return sum;
    }

    @Override
    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        if (!(potential instanceof P1ImageHarmonic)) return;
        sum += ((P1ImageHarmonic)potential).getDUDW(atoms);
    }
}
