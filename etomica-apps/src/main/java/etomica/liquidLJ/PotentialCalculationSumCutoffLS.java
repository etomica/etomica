/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.liquidLJ;

import etomica.atom.IAtomList;
import etomica.potential.IPotentialAtomic;
import etomica.potential.PotentialCalculation;

/**
 * Evaluates the energy summed over all iterated atoms. Each call to doCalculate
 * accumulates the energy of the given potential applied to the atoms produced
 * by the given iterator. Each such call accumulates a single sum, which can
 * be accessed via the getSum() method.  Sum is re-zeroed only upon a call to the
 * zeroSum() method. 
 *
 * @author David Kofke
 */
public class PotentialCalculationSumCutoffLS implements PotentialCalculation {

    public PotentialCalculationSumCutoffLS() {
        sums = new double[2][0];
    }
    
    /**
	 * Adds to the energy sum the energy values obtained from application of the given potential to the
	 * atoms produced by the given iterator.  Iterator is reset by method before beginning calculation.
	 */
	public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
	    if (!(potential instanceof Potential2SoftSphericalLSMulti)) return;
	    Potential2SoftSphericalLSMulti p2 = (Potential2SoftSphericalLSMulti)potential;
        double[][] ijSums = p2.energyVirialCut(atoms);
        int n = ijSums[0].length;
        if (n!= sums[0].length) {
            sums[0] = new double[n];
            sums[1] = new double[n];
        }
        for (int i=0; i<sums[0].length; i++) {
            sums[0][i] += ijSums[0][i];
            sums[1][i] += ijSums[1][i];
        }
	}
	
	/**
	 * Sets the energy sum to zero, typically to begin a new energy-sum calculation.
	 */
	public void zeroSums() {
	    for (int i=0; i<sums[0].length; i++) {
	        sums[0][i] = sums[1][i] = 0.0;
	    }
	}

	/**
	 * Returns the current value of the energy sum.
	 */
	public double[][] getSums() {
        return sums;
    }
	
	protected final double[][] sums;
}
