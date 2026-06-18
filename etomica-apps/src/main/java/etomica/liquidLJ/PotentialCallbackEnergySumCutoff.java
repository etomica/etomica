/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.liquidLJ;

import etomica.potential.compute.PotentialCallback;
import etomica.space.Vector;

import java.util.Arrays;

/**
 * Evaluates the energy summed over all iterated atoms. Each call to doCalculate
 * accumulates the energy of the given potential applied to the atoms produced
 * by the given iterator. Each such call accumulates a single sum, which can
 * be accessed via the getSum() method.  Sum is re-zeroed only upon a call to the
 * zeroSum() method. 
 *
 * @author David Kofke
 */
public class PotentialCallbackEnergySumCutoff implements PotentialCallback {

    public PotentialCallbackEnergySumCutoff(double[] cutoffs) {
        r2Cuts = new double[cutoffs.length];
        for (int i=0; i<cutoffs.length; i++) {
            r2Cuts[i] = cutoffs[i]*cutoffs[i];
        }
        sums = new double[cutoffs.length];
    }
    
    /**
	 * Adds to the energy sum the energy values obtained from application of the given potential to the
	 * atoms produced by the given iterator.  Iterator is reset by method before beginning calculation.
	 */
    @Override
    public void pairCompute(int iAtom, int jAtom, Vector dr, double[] u012) {
        double r2 = dr.squared();
        for (int i=sums.length-1; i>=0; i--) {
            if (r2 > r2Cuts[i]) break;
            sums[i] += u012[0];
        }
	}
	
	/**
	 * Sets the energy sum to zero, typically to begin a new energy-sum calculation.
	 */
	public void zeroSum() {
        Arrays.fill(sums, 0);
	}

	/**
	 * Returns the current value of the energy sum.
	 */
	public double[] getSums() {
        return sums;
    }
	
	protected double[] sums, r2Cuts;
}
