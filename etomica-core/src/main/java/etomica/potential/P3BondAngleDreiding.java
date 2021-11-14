/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Simple 3-body soft bond-angle for Dreiding potential 
 * 
 *  * U(theta) = gamma * [ cos(theta) - cos(theta_e)]^2
 * 
 *  
 * where gamma is the potential's energy parameter
 * 		 theta_e is the equilibrium bond angle
 * 
 * 
 * @author Tai Tan
 */

public class P3BondAngleDreiding implements PotentialSoft {
	
	public P3BondAngleDreiding (Space space, double thetaEq, double gamma){
		dr12 = space.makeVector();
		dr23 = space.makeVector();
		setAngle(thetaEq);
		setGamma(gamma);
		gradient = new Vector[3];
		gradient[0] = space.makeVector();
		gradient[1] = space.makeVector();
		gradient[2] = space.makeVector();
	}

	public double energy(IAtomList atomSet){
        IAtom atom0 = atomSet.get(0);
        IAtom atom1 = atomSet.get(1);
        IAtom atom2 = atomSet.get(2);
		dr12.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
		dr23.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
		boundary.nearestImage(dr12);
		boundary.nearestImage(dr23);
		
		double costheta = -dr12.dot(dr23)/Math.sqrt(dr12.squared()*dr23.squared());
		
		//machine precision can give us numbers with magnitudes slightly greater than 1
		if (costheta > 1){
			costheta = 1;
		}
		
        else if (costheta < -1) {
            costheta = -1;
        }
		
        double dtheta = costheta - Math.cos(thetaEq); 
        return gamma*dtheta*dtheta;
		
	}

	 /**
     * Sets the nominal bond angle (in radians)
     */
    public void setAngle(double newAngle) {
        thetaEq = newAngle;
    }
    
    /**
     * Returns the nominal bond angle (in radians)
     */
    public double getAngle() {
        return thetaEq;
    }

    /**
     * Sets the characteristic energy of the potential
     */
    public void setGamma(double newGamma) {
        gamma = newGamma;
    }
    
    /**
     * Returns the characteristic energy of the potential
     */
    public double getGamma() {
        return gamma;
    }
    
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }
	
    protected final Vector dr12, dr23;
    protected final Vector[] gradient;
    protected Boundary boundary;
    private double gamma;
    private double thetaEq;
    private static final long serialVersionUID = 1L;
}
