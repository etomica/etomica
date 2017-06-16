/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.space.Space;

/**
 * Simple electrostatic potential class.
 * @author Andrew Schultz
 */
public class P2ElectrostaticWithHardCore extends Potential2SoftSpherical {

    public P2ElectrostaticWithHardCore(Space space) {
        super(space);
    }
    
    public double d2u(double r2) {
        return +2*u(r2);
    }

    public double du(double r2) {
        return -u(r2);
    }

    public double uInt(double rc) {
        // lie.  Nobody really wants to know it's infinity
        return 0;
    }

    public double u(double r2) {
    	double r = Math.sqrt(r2);
    	/*double sigma;
    	if (qq > 0) {
    		sigma = 0;
    	} else {
    		sigma = this.sigma;
    	}
    	*/
    	
    	
    	if (r > sigma) {
    		/*if (checkE) {
    			System.out.println("electrostatic potential is " + qq/r + " for r of " + r + ", q1 of " + q1 + ", and q2 of " + q2);
    		}*/
    		return qq/r;
    		
    	} else {
    		//System.out.println("electrostatic potential is inf");
    		return Double.POSITIVE_INFINITY;
    		//return 0;
    	}
    }

    /**
     * Returns the charge on the first type of atom.
     */
    public double getCharge1() {
        return q1;
    }
    
    /**
     * Sets the charge on the first type of atom.
     */
    public void setCharge1(double newCharge1) {
        q1 = newCharge1;
        qq = q1*q2;
    }

    /**
     * Returns the charge on the second type of atom.
     */
    public double getCharge2() {
        return q2;
    }
    
    /**
     * Sets the charge on the second type of atom.
     */
    public void setCharge2(double newCharge2) {
        q2 = newCharge2;
        qq = q1*q2;
    }
    
    public void setSigma(double sigma) {
    	this.sigma = sigma;
    }
    
    private static final long serialVersionUID = 1L;
    protected double q1, q2, qq;
    protected double sigma = 0;
    public static boolean checkE;
}

