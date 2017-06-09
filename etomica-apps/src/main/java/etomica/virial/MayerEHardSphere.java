/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.box.Box;
import etomica.atom.IMoleculeList;

/**
 * @author kofke
 *
 * The e-function for hard spheres, returning 0 for r<sigma, 1 otherwise.
 */
public class MayerEHardSphere extends MayerESpherical {

	private double sigma, sigma2;

	/**
	 * Constructor for MayerEHardSphere.
	 * @param potential
	 */
	public MayerEHardSphere() {
		this(1.0);
	}
	
	public MayerEHardSphere(double sigma) {
		super(null);
		setSigma(sigma);
	}

	public double f(IMoleculeList pair, double r2, double beta) {
		return (r2<sigma2) ? 0.0 : 1.0;
	}

	/**
	 * Returns the HS diameter.
	 * @return double
	 */
	public double getSigma() {
		return sigma;
	}

	/**
	 * Sets the HS diameter.
	 * @param sigma The sigma to set
	 */
	public void setSigma(double sigma) {
		this.sigma = sigma;
		sigma2 = sigma*sigma;
	}

    public void setBox(Box newBox) {
    }
}
