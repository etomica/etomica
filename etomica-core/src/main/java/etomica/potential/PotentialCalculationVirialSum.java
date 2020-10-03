/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomList;

/**
 * Evaluates the virial summed over all iterated atoms.
 *
 * @author David Kofke
 */
public class PotentialCalculationVirialSum implements PotentialCalculation {
		
    /**
	 * Adds to the virial sum the energy values obtained from application of the given potential to the
	 * atoms produced by the given iterator.  Iterator is reset by method before beginning calculation.
	 */
	public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        if (!(potential instanceof PotentialSoft)) {
            return;
        }
        sum += ((PotentialSoft)potential).virial(atoms);	// This will sum the energy values.
	}
	
	/**
	 * Sets the virial sum to zero, typically to begin a new virial-sum calculation.
	 * @return this instance, so the method can be called in-line as the instance is
	 * passed to the PotentialMaster.
	 */
	public PotentialCalculationVirialSum zeroSum() {
		sum = 0.0;
		return this;
	}

    private static final long serialVersionUID = 1L;

	/**
	 * Returns the current value of the energy sum.
	 */
	public double getSum() {return sum;}
	
	private double sum = 0.0;

 }//end VirialSum
