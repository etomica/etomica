/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtomList;
import etomica.normalmode.Potential2SoftSphericalLSMultiLat.ReturnValue;
import etomica.potential.IPotentialAtomic;
import etomica.space.Space;

/**
 * Sums the force on each iterated atom and adds it to the integrator agent
 * associated with the atom.
 */
public class PotentialCalculationSolidSuperCutLS extends PotentialCalculationSolidSuperCut {
        
    public PotentialCalculationSolidSuperCutLS(Space space, CoordinateDefinition coordinateDefinition, double[] cutoffs) {
        super(space, coordinateDefinition, cutoffs);
    }
    
    /**
     * Adds forces due to given potential acting on the atoms produced by the iterator.
     * Implemented for only 1- and 2-body potentials.
     */
    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        if (!(potential instanceof Potential2SoftSphericalLSMultiLat)) return;
        Potential2SoftSphericalLSMultiLat potentialSoft = (Potential2SoftSphericalLSMultiLat)potential;
        ReturnValue rv = potentialSoft.energyVirialCut(atoms);
        int n = rv.energySum.length;
        if (n != energySum.length) {
            energySum = new double[n];
            virialSum = new double[n];
            sum1 = new double[n];
            dadbSum = new double[n];
            pzxySum = new double[n];
        }
        for (int i=0; i<n; i++) {
            energySum[i] += rv.energySum[i];
            virialSum[i] += rv.virialSum[i];
            sum1[i] += rv.sum1[i]*fac1;
            dadbSum[i] += rv.dadbSum[i];
            pzxySum[i] += rv.pzxySum[i];
        }
    }
}
