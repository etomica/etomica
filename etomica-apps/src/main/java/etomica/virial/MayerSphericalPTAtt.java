/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.box.Box;
import etomica.math.SpecialFunctions;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotential;
import etomica.potential.Potential2Spherical;

/**
 * Mayer function class that returns temperature expansion bond of the desired
 * order.
 * 
 * @author Andrew Schultz
 */
public class MayerSphericalPTAtt implements MayerFunction {

	/**
	 * Constructor takes reference potential, attractive component and order in temperature.
	 */
	public MayerSphericalPTAtt(Potential2Spherical potentialRef, Potential2Spherical potentialAtt, int order) {
		this.potentialRef = potentialRef;
		this.potentialAtt = potentialAtt;
		this.order = order;
		this.factOrder = SpecialFunctions.factorial(order);
	}

	public double f(IMoleculeList pair, double r2, double beta) {
		double uRef = potentialRef.u(r2);
		if (Double.isInfinite(uRef)) {
			return 0;
		}
		// double dfdkT = Math.exp(-beta*uRef)*Math.pow(-potentialAtt.u(r2)*beta, order)/SpecialFunctions.factorial(order);
		double dfdkT = Math.exp(-beta*uRef)/factOrder;
		if (order > 0) {
		    double betaUA = -potentialAtt.u(r2)*beta;
		    switch (order) {
		        case 1: 
		            dfdkT *= betaUA;
		            break;
		        case 2: 
		            dfdkT *= betaUA*betaUA;
		            break;
                case 3: 
                    dfdkT *= betaUA*betaUA*betaUA;
                    break;
                case 4:
                    double betaUA2 = betaUA*betaUA;
                    dfdkT *= betaUA2*betaUA2;
                    break;
                case 5:
                    betaUA2 = betaUA*betaUA;
                    dfdkT *= betaUA2*betaUA2*betaUA;
                    break;
                case 6:
                    double betaUA3 = betaUA*betaUA*betaUA;
                    dfdkT *= betaUA3*betaUA3;
                    break;
                case 7:
                    betaUA3 = betaUA*betaUA*betaUA;
                    dfdkT *= betaUA3*betaUA3*betaUA;
                    break;
                case 8:
                    betaUA2 = betaUA*betaUA;
                    double betaUA4 = betaUA2*betaUA2;
                    dfdkT *= betaUA4*betaUA4;
                    break;
                case 9:
                    betaUA3 = betaUA*betaUA*betaUA;
                    dfdkT *= betaUA3*betaUA3*betaUA3;
                    break;
                default:
                    throw new RuntimeException("oops, please implement order "+order);
		    }
		}
		
		
		if (Double.isNaN(dfdkT)) {
			throw new RuntimeException ("dfdT is NaN");
		}
		return dfdkT;
	}
	
	public void setBox(Box newBox) {
	    potentialRef.setBox(newBox);
	    potentialAtt.setBox(newBox);
	}

    public IPotential getPotential() {
        return potentialRef;
    }

    protected final Potential2Spherical potentialRef, potentialAtt;
	protected final int order;
	protected final long factOrder;
}
