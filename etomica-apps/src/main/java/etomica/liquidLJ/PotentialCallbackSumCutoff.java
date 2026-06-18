/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.liquidLJ;

import etomica.potential.compute.PotentialCallback;
import etomica.space.Vector;

import java.util.Arrays;

public class PotentialCallbackSumCutoff implements PotentialCallback {

    public PotentialCallbackSumCutoff(double[] cutoffs) {
        r2Cuts = new double[cutoffs.length];
        for (int i=0; i<cutoffs.length; i++) {
            r2Cuts[i] = cutoffs[i]*cutoffs[i];
        }
        uSums = new double[cutoffs.length];
        vSums = new double[cutoffs.length];
    }

    /**
	 * Adds to the energy sum the energy values obtained from application of the given potential to the
	 * atoms produced by the given iterator.  Iterator is reset by method before beginning calculation.
	 */
    @Override
    public void pairCompute(int iAtom, int jAtom, Vector dr, double[] u012) {
        double r2 = dr.squared();
        if (r2 > r2Cuts[r2Cuts.length-1]) return;
        for (int i=uSums.length-1; i>=0; i--) {
            if (r2 > r2Cuts[i]) break;
            uSums[i] += u012[0];
            vSums[i] += u012[1];
        }
	}
	
	/**
	 * Sets the energy sum to zero, typically to begin a new energy-sum calculation.
	 */
	public void zeroSums() {
	    Arrays.fill(uSums, 0);
        Arrays.fill(vSums, 0);
	}

	/**
	 * Returns the current value of the energy sum.
	 */
	public double[] getUSums() {
        return uSums;
    }
	
	/**
     * Returns the current value of the energy sum.
     */
    public double[] getVSums() {
        return vSums;
    }

	
	protected double[] uSums, vSums, r2Cuts;
}
