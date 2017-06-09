/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.box.Box;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotential;

/**
 * @author kofke
 *
 * Hard-sphere Mayer function.  -1 if r < sigma; 0 otherwise
 */
public class MayerHardSphere implements MayerFunction {

    private static final long serialVersionUID = 1L;
    private double sigma, sigma2;
    /**
     * Constructor for MayerHardSphere.
     */
    public MayerHardSphere() {
        this(1.0);
    }
    public MayerHardSphere(double sigma) {
        setSigma(sigma);
    }

    /**
     * @see etomica.virial.MayerFunctionSpherical#f(etomica.AtomPair)
     */
    public double f(IMoleculeList pair, double r2, double beta) {
        return (r2<sigma2) ? -1.0 : 0.0;
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

    public IPotential getPotential() {
        return null;
    }

    public void setBox(Box newBox) {
    }
}
