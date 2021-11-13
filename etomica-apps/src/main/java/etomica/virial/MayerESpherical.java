/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.box.Box;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotential;
import etomica.potential.Potential2Spherical;
/**
 * @author kofke
 *
 * Simple e-bond, exp(-beta*u).
 */
public class MayerESpherical implements MayerFunction {

	/**
	 * Constructor for MayerESpherical.
	 */
	public MayerESpherical(Potential2Spherical potential) {
		this.potential = potential;
	}

	/**
	 * @see etomica.virial.MayerFunctionSpherical#f(etomica.AtomPair, double, double)
	 */
	public double f(IMoleculeList pair, double r2, double beta) {
		return Math.exp(-beta*potential.u(r2));
	}
	
	public void setBox(Box newBox) {
	}

	private final Potential2Spherical potential;

	public IPotential getPotential() {
		return potential;
	}

}
