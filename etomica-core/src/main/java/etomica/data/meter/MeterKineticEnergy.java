/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.units.dimensions.Energy;

/**
 * Meter for the total kinetic energy in a box Computes total KE by summing values of KE returned by every atom in the
 * box. A different box-dependent atom integrator may be set to permit calculation over a particular set of atoms in the
 * box.
 */
public class MeterKineticEnergy extends DataSourceScalar {

    private final Box box;

    public MeterKineticEnergy(Box box) {
        super("Kinetic Energy", Energy.DIMENSION);
        this.box = box;
    }

    /**
     * Returns the total kinetic energy summed over all atoms produced by the iterator when applied to the given box.
     * Does not include contributions from atoms having infinite mass (it assumes they are stationary).
     */
    public double getDataAsScalar() {
        double ke = 0.0;
        for (IAtom atom : box.getLeafList()) {
            double mass = atom.getType().getMass();
            if (mass == Double.POSITIVE_INFINITY) continue;
            ke += 0.5 * mass * (((IAtomKinetic)atom).getVelocity().squared());
        }
        return ke;
    }

    /**
     * @return Returns the box.
     */
    public Box getBox() {
        return box;
    }
}
