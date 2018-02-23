/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.liquidLJ;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.IPotentialAtomic;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialCalculation;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Evaluates the energy summed over all iterated atoms. Each call to doCalculate
 * accumulates the energy of the given potential applied to the atoms produced
 * by the given iterator. Each such call accumulates a single sum, which can
 * be accessed via the getSum() method.  Sum is re-zeroed only upon a call to the
 * zeroSum() method. 
 *
 * @author David Kofke
 */
public class PotentialCalculationEnergySumCutoff implements PotentialCalculation {

    public PotentialCalculationEnergySumCutoff(Space space, double[] cutoffs) {
        dr = space.makeVector();
        r2Cuts = new double[cutoffs.length];
        for (int i=0; i<cutoffs.length; i++) {
            r2Cuts[i] = cutoffs[i]*cutoffs[i];
        }
        sums = new double[cutoffs.length];
    }
    
    public void setBox(Box box) {
        this.box = box;
        boundary = box.getBoundary();
    }

    /**
	 * Adds to the energy sum the energy values obtained from application of the given potential to the
	 * atoms produced by the given iterator.  Iterator is reset by method before beginning calculation.
	 */
	public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        dr.Ev1Mv2(atoms.get(1).getPosition(),atoms.get(0).getPosition());
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        double u = ((Potential2SoftSpherical)potential).u(r2);
        for (int i=sums.length-1; i>=0; i--) {
            if (r2 > r2Cuts[i]) break;
            sums[i] += u;
        }
	}
	
	/**
	 * Sets the energy sum to zero, typically to begin a new energy-sum calculation.
	 */
	public void zeroSum() {
	    for (int i=0; i<sums.length; i++) {
	        sums[i] = 0.0;
	    }
		boundary = box.getBoundary();
	}

	/**
	 * Returns the current value of the energy sum.
	 */
	public double[] getSums() {
        return sums;
    }
	
	protected double[] sums, r2Cuts;
	protected final Vector dr;
	protected Box box;
	protected Boundary boundary;
}
