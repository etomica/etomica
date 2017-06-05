/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomList;
import etomica.api.IMoleculeList;
import etomica.api.IPotentialAtomic;
import etomica.api.IPotentialMolecular;

/**
 * Evaluates the virial summed over all iterated molecules.
 *
 * @author Tai Boon Tan
 */
public class PotentialCalculationMolecularVirialSum implements PotentialCalculationMolecular {
		
	
	public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
		// TODO Auto-generated method stub
		
	}
	
	public void doCalculation(IMoleculeList molecules, IPotentialMolecular potential) {
		if (!(potential instanceof PotentialMolecularSoft)) {
			return;
		}	
		sum += ((PotentialMolecularSoft)potential).virial(molecules);
	}
	
	/**
	 * Sets the virial sum to zero, typically to begin a new virial-sum calculation.
	 * @return this instance, so the method can be called in-line as the instance is
	 * passed to the PotentialMaster.
	 */
	public PotentialCalculationMolecularVirialSum zeroSum() {
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
