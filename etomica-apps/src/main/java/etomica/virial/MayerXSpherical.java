/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.api.IPotential;
import etomica.box.Box;
import etomica.molecule.IMoleculeList;
import etomica.potential.Potential2SoftSpherical;
/**
 * @author kate
 *
 * -beta*e12*r12*du12/dr12 bond
 */
public class MayerXSpherical implements MayerFunction {

	/**
	 * Constructor for MayerESpherical.
	 */
	public MayerXSpherical(Potential2SoftSpherical potential) {
		this.potential = potential;
	}

	/**
	 * @see etomica.virial.MayerFunctionSpherical#f(etomica.AtomPair, double, double)
	 */
	public double f(IMoleculeList pair, double r2, double beta) {

		if (r2 == 0) {
			return 0;
		} else {
			return Math.exp(-beta*potential.u(r2))*(-beta*potential.du(r2));
		}
	}
	
	public void setBox(Box newBox) {
	    potential.setBox(newBox);
	}

	private final Potential2SoftSpherical potential;

	public IPotential getPotential() {
		return potential;
	}

}
