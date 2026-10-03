/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;

/**
 * 3-body soft bond-angle potential with tunable stiffness
 *
 * @author andrew
 */
public class P3BondAngleStiffChain implements IPotentialBondAngle {

    public P3BondAngleStiffChain(double epsilon) {
        setEpsilon(epsilon);
    }

    public double u(double costheta) {
        return epsilon*(1+costheta);
    }

    /**
     * Sets the characteristic energy of the potential
     */
    public void setEpsilon(double newEpsilon) {
        epsilon = newEpsilon;
    }
    
    /**
     * Returns the characteristic energy of the potential
     */
    public double getEpsilon() {
        return epsilon;
    }
    
    public Dimension getEpsilonDimension() {
        return Energy.DIMENSION;
    }
    
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public void udu(double costheta, double[] u, double[] du) {
        du[0] = epsilon;
        u[0] = epsilon*(1 + costheta);
    }

    protected double epsilon;
}
