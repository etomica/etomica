/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.meam;

import etomica.atom.IAtomList;
import etomica.potential.IPotentialAtomic;
import etomica.potential.PotentialCalculationEnergySum;

/**
 * Evaluate the energy for the EAM model using the P2EAM potential.  This
 * potential calculation will explicitly include the embedding contribution
 * from the EAM potential.
 */
public class PotentialCalculationEnergySumEAM extends PotentialCalculationEnergySum {
    
    protected final P2EAM p2;
    protected double sum;
    
    public PotentialCalculationEnergySumEAM(P2EAM p2) {
        this.p2 = p2;
    }
    
    /**
     * Adds to the energy sum the energy values obtained from application of the given potential to the
     * atoms.
     */
    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        sum += potential.energy(atoms);
    }
    
    /**
     * Sets the energy sum to zero, typically to begin a new energy-sum calculation.
     */
    public void zeroSum() {
        p2.reset();
        sum = 0.0;
    }

    /**
     * Returns the current value of the energy sum.
     */
    public double getSum() {
        return sum + p2.energy1();
    }
}
