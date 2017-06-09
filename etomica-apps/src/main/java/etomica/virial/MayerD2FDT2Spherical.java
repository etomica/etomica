/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.box.Box;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotential;
import etomica.potential.Potential2Spherical;
/**
 * Required for computing second derivatives  of virial coefficients w/r to temperature
 * @author kate
 */
public class MayerD2FDT2Spherical implements MayerFunction {

	/**
	 * Constructor for MayerESpherical.
	 */
	public MayerD2FDT2Spherical(Potential2Spherical potential) {
		this.potential = potential;
	}

	/**
	 * @see etomica.virial.MayerFunctionSpherical#f(etomica.AtomPair, double, double)
	 */
	public double f(IMoleculeList pair, double r2, double beta) {
		double u = potential.u(r2);
		if (Double.isInfinite(u)) {
			return 0;
		}
		
		
		double dfdkT = Math.exp(-beta*u)*u*beta*beta;
		
		double d2fdkT2 = dfdkT*(-2.0*beta + u*beta*beta);
		
		
		return d2fdkT2;
	}
	
	public void setBox(Box newBox) {
	    potential.setBox(newBox);
	}

	private final Potential2Spherical potential;

	public IPotential getPotential() {
		return potential;
	}

}
