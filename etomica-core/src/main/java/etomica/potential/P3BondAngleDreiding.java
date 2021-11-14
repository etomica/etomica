/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Tensor;
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

	public Vector[] gradient(IAtomList atomSet, Tensor pressureTensor){
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
       //double diff = Math.acos(costheta)-thetaEq;
       //System.out.println(Math.acos(costheta)+" "+ thetaEq + " "+ diff/thetaEq*100);
		
		double constant = 2*gamma*dtheta;
		double dr12sq = Math.sqrt(dr12.squared()); 
		double dr23sq = Math.sqrt(dr23.squared());
		double dotProd = dr12.dot(dr23);
		
		
		/*
		 *  Atom 1
		 */
		double numdx1a   =  (dr23.getX(0));
		double denomdx1a = dr12sq * dr23sq;
		double numdx1b   = -(dr12.getX(0))* dotProd;
		double denomdx1b = dr12.squared()*dr12sq*dr23sq;
		gradient[0].setX(0, numdx1a / denomdx1a + numdx1b / denomdx1b);
		
		double numdy1a   =  (dr23.getX(1));
		double denomdy1a = dr12sq * dr23sq;
		double numdy1b   = -(dr12.getX(1))* dotProd;
		double denomdy1b = dr12.squared()*dr12sq*dr23sq;
		gradient[0].setX(1, numdy1a / denomdy1a + numdy1b / denomdy1b);
		
		double numdz1a   =  (dr23.getX(2));
		double denomdz1a = dr12sq * dr23sq;
		double numdz1b   = -(dr12.getX(2))* dotProd;
		double denomdz1b = dr12.squared()*dr12sq*dr23sq;
		gradient[0].setX(2, numdz1a / denomdz1a + numdz1b / denomdz1b);
		
		
		/*
		 *  Atom 2
		 */
		double numdx2a   =   ((dr12.getX(0)) - (dr23.getX(0)));
		double denomdx2a = dr12sq * dr23sq;
		double numdx2b   =   (dr12.getX(0)) * dotProd;
		double denomdx2b = dr12.squared()*dr12sq*dr23sq;
		double numdx2c   =  -(dr23.getX(0)) * dotProd;
		double denomdx2c = dr12sq*dr23sq*dr23.squared();
		gradient[1].setX(0, numdx2a / denomdx2a + numdx2b / denomdx2b + numdx2c / denomdx2c);
		
		double numdy2a   =  ((dr12.getX(1)) - (dr23.getX(1)));
		double denomdy2a = dr12sq * dr23sq;
		double numdy2b   =  (dr12.getX(1)) * dotProd;
		double denomdy2b = dr12.squared()*dr12sq*dr23sq;
		double numdy2c   =  -(dr23.getX(1)) * dotProd;
		double denomdy2c = dr12sq*dr23sq*dr23.squared();
		gradient[1].setX(1, numdy2a / denomdy2a + numdy2b / denomdy2b + numdy2c / denomdy2c);
		
		double numdz2a   =  ((dr12.getX(2)) - (dr23.getX(2)));
		double denomdz2a = dr12sq * dr23sq;
		double numdz2b   =  (dr12.getX(2)) * dotProd;
		double denomdz2b = dr12.squared()*dr12sq*dr23sq;
		double numdz2c   =  -(dr23.getX(2)) * dotProd;
		double denomdz2c = dr12sq*dr23sq*dr23.squared();
		gradient[1].setX(2, numdz2a / denomdz2a + numdz2b / denomdz2b + numdz2c / denomdz2c);
		
		
        /*
         *  Atom 3
         */
        double numdx3a   =  -(dr12.getX(0));  
        double denomdx3a = dr12sq*dr23sq;
        double numdx3b   =   (dr23.getX(0))*dotProd;
        double denomdx3b = dr12sq*dr23sq*dr23.squared();
        gradient[2].setX(0, numdx3a / denomdx3a + numdx3b / denomdx3b);
        
        double numdy3a   =  -(dr12.getX(1));  
        double denomdy3a = dr12sq*dr23sq;
        double numdy3b   =   (dr23.getX(1))*dotProd;
        double denomdy3b = dr12sq*dr23sq*dr23.squared();
        gradient[2].setX(1, numdy3a / denomdy3a + numdy3b / denomdy3b);
        
        double numdz3a   =  -(dr12.getX(2));  
        double denomdz3a = dr12sq*dr23sq;
        double numdz3b   =   (dr23.getX(2))*dotProd;
        double denomdz3b = dr12sq*dr23sq*dr23.squared();
        gradient[2].setX(2, numdz3a / denomdz3a + numdz3b / denomdz3b);
        
        gradient[0].TE(constant);
        gradient[1].TE(constant);
        gradient[2].TE(constant);
        
        if (pressureTensor != null){
   //     	pressureTensor.PEv1v2(gradient[0],dr);
        }
		return gradient;
	}
	
    public Vector[] gradient(IAtomList atoms) {

        return gradient(atoms,null);
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
