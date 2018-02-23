/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Boundary;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space3d.Vector3D;


/**
 * Simple 4-body soft bond-angle for Dreiding potential 
 * 
 *   U(phi) = (1/2) * Vjk * {1 - cos [njk (phi - phiEq)]}
 *   
 * where Vjk is the rotation barrier (always +ve) [unit Kelvin]
 * 		 njk is the periodicity (an integer)
 * 		 phiEq is the equilibrium dihedral angle  [unit radians]
 * 
 * @author Tai Tan
 */

public class P4TorsionDreiding extends Potential implements PotentialSoft {
	
	public P4TorsionDreiding (Space space, double phiEq, double beta, int n){
		super(4, space);
		
		dr12 = new Vector3D();
		dr23 = new Vector3D();
		dr34 = new Vector3D();
		dra  = new Vector3D();
		drb  = new Vector3D();
		drbxa = new Vector3D();
		
		setPhi(phiEq);
		setBeta(beta);
		setN(n);
		
		gradient = new Vector[4];
		gradient[0] = space.makeVector();
		gradient[1] = space.makeVector();
		gradient[2] = space.makeVector();
		gradient[3] = space.makeVector();
		}
	
	public double virial(IAtomList atomSet){
		return 0.0;
	}
	
	public Vector[] gradient(IAtomList atomSet, Tensor pressureTensor){
        IAtom atom0 = atomSet.get(0);
        IAtom atom1 = atomSet.get(1);
        IAtom atom2 = atomSet.get(2);
        IAtom atom3 = atomSet.get(3);
		dr12.Ev1Mv2(atom0.getPosition(), atom1.getPosition());
		dr23.Ev1Mv2(atom1.getPosition(), atom2.getPosition());
		dr34.Ev1Mv2(atom3.getPosition(), atom2.getPosition());
		boundary.nearestImage(dr12);
		boundary.nearestImage(dr23);
		boundary.nearestImage(dr34);
		
		dra.E(dr12);
		drb.E(dr34);
		
		dra.XE(dr23);
		drb.XE(dr23);
		
		drbxa.E(drb);
		drbxa.XE(dra);

		/*
		 *  Atom 1
		 */
		double const1 = -Math.sqrt(dr23.squared()) / dra.squared();
		gradient[0].setX(0, const1 * dra.getX(0));
		gradient[0].setX(1, const1 * dra.getX(1));
		gradient[0].setX(2, const1 * dra.getX(2));
		
		/*
		 *  Atom 2
		 */
		double const2a = Math.sqrt(dr23.squared()) / dra.squared();
		double const2b = dr12.dot(dr23) / (dra.squared()*Math.sqrt(dr23.squared()));
		double const2c = dr34.dot(dr23) / (drb.squared()*Math.sqrt(dr23.squared()));
		gradient[1].setX(0, const2a*dra.getX(0) + const2b*dra.getX(0) - const2c*drb.getX(0));
		gradient[1].setX(1, const2a*dra.getX(1) + const2b*dra.getX(1) - const2c*drb.getX(1));
		gradient[1].setX(2, const2a*dra.getX(2) + const2b*dra.getX(2) - const2c*drb.getX(2));
		
        /*
         *  Atom 3
         */
		double const3a = Math.sqrt(dr23.squared()) / drb.squared();
		double const3b = dr34.dot(dr23) / (drb.squared()*Math.sqrt(dr23.squared()));
		double const3c = dr12.dot(dr23) / (dra.squared()*Math.sqrt(dr23.squared()));
		gradient[2].setX(0, - const3a*drb.getX(0) + const3b*drb.getX(0) - const3c*dra.getX(0));
		gradient[2].setX(1, - const3a*drb.getX(1) + const3b*drb.getX(1) - const3c*dra.getX(1));
		gradient[2].setX(2, - const3a*drb.getX(2) + const3b*drb.getX(2) - const3c*dra.getX(2));
		
		/*
		 *  Atom 4
		 */
		double const4 = -Math.sqrt(dr23.squared()) / drb.squared();
		gradient[3].setX(0, const4 * drb.getX(0));
		gradient[3].setX(1, const4 * drb.getX(1));
		gradient[3].setX(2, const4 * drb.getX(2));
		
		return gradient;
	}
	
    public Vector[] gradient(IAtomList atoms) {

        return gradient(atoms,null);
    }
	
	
	public void setBox(Box box){
		boundary = box.getBoundary();
	}
	
	public double energy(IAtomList atomSet){
        IAtom atom0 = atomSet.get(0);
        IAtom atom1 = atomSet.get(1);
        IAtom atom2 = atomSet.get(2);
        IAtom atom3 = atomSet.get(3);
		dr12.Ev1Mv2(atom0.getPosition(), atom1.getPosition());
		dra. Ev1Mv2(atom0.getPosition(), atom1.getPosition());
		dr23.Ev1Mv2(atom1.getPosition(), atom2.getPosition());
		dr34.Ev1Mv2(atom3.getPosition(), atom2.getPosition());
		drb. Ev1Mv2(atom0.getPosition(), atom1.getPosition());
		
		boundary.nearestImage(dr12);
		boundary.nearestImage(dr23);
		boundary.nearestImage(dr34);
		
		/*
		 * To get the torsional angle
		 */
		dra.XE(dr23);
		drb.XE(dr23);

		double cosphi = dra.dot(drb)/Math.sqrt(dra.squared()*drb.squared());
		double phi;
		
		//machine precision can give us numbers with magnitudes slightly greater than 1
		if (cosphi > 1){
			phi = 0;
		}
		
        else if (cosphi < -1) {
            phi = Math.PI;
        }
        else {
            phi = Math.acos(cosphi);
        }
		 
        return beta*(1-Math.cos(n*(phi - phiEq)));
	}
	
	/*
	 *  First Derivative of energy du         
	 */
	
	public double du(IAtomList atomSet){
        IAtom atom0 = atomSet.get(0);
        IAtom atom1 = atomSet.get(1);
        IAtom atom2 = atomSet.get(2);
        IAtom atom3 = atomSet.get(3);
		dr12.Ev1Mv2(atom0.getPosition(), atom1.getPosition());
		dr23.Ev1Mv2(atom1.getPosition(), atom2.getPosition());
		dr34.Ev1Mv2(atom3.getPosition(), atom2.getPosition());
		boundary.nearestImage(dr12);
		boundary.nearestImage(dr23);
		boundary.nearestImage(dr34);
		
		dra.E(dr12);
		drb.E(dr34);
		
		dra.XE(dr23);
		drb.XE(dr23);
		
		drbxa.E(drb);
		drbxa.XE(dra);

		double cosphi = dra.dot(drb)/Math.sqrt(dra.squared()*drb.squared());
		double sinphi = (drbxa.dot(dr23))/Math.sqrt(dra.squared()*drb.squared()*dr23.squared());
				
		double dx1, dy1, dz1;
		double dx2, dy2, dz2;
		double dx3, dy3, dz3;
		double dx4, dy4, dz4;
		double du1, du2, du3, du4;
		
		/*
		 *  Atom 1
		 */
		double const1 = -Math.sqrt(dr23.squared()) / dra.squared();
		dx1 = const1 * dra.getX(0);
		dy1 = const1 * dra.getX(1);
		dz1 = const1 * dra.getX(2);
		du1 = dx1 + dy1 + dz1;
		
		/*
		 *  Atom 2
		 */
		double const2a = Math.sqrt(dr23.squared()) / dra.squared();
		double const2b = dr12.dot(dr23) / (dra.squared()*Math.sqrt(dr23.squared()));
		double const2c = dr34.dot(dr23) / (drb.squared()*Math.sqrt(dr23.squared()));
		dx2 = const2a*dra.getX(0) + const2b*dra.getX(0) - const2c*drb.getX(0);
		dy2 = const2a*dra.getX(1) + const2b*dra.getX(1) - const2c*drb.getX(1);
		dz2 = const2a*dra.getX(2) + const2b*dra.getX(2) - const2c*drb.getX(2);
		du2 = dx2 + dy2 + dz2;
		
        /*
         *  Atom 3
         */
		double const3a = Math.sqrt(dr23.squared()) / drb.squared();
		double const3b = dr34.dot(dr23) / (drb.squared()*Math.sqrt(dr23.squared()));
		double const3c = dr12.dot(dr23) / (dra.squared()*Math.sqrt(dr23.squared()));
		dx3 = - const3a*drb.getX(0) + const3b*drb.getX(0) - const3c*dra.getX(0);
		dy3 = - const3a*drb.getX(1) + const3b*drb.getX(1) - const3c*dra.getX(1);
		dz3 = - const3a*drb.getX(2) + const3b*drb.getX(2) - const3c*dra.getX(2);
		du3 = dx3 + dy3 + dz3;
		
		/*
		 *  Atom 4
		 */
		double const4 = -Math.sqrt(dr23.squared()) / drb.squared();
		dx4 = const4 * drb.getX(0);
		dy4 = const4 * drb.getX(1);
		dz4 = const4 * drb.getX(2);
		du4 = dx4 + dy4 + dz4;
		
		return dudphi(n, cosphi, sinphi)* (du1 + du2 + du3 + du4);
	}

	
    public void setPhi(double newPhi) {
        phiEq = newPhi;
    }

    public double getPhi() {
        return phiEq;
    }

    public void setBeta(double newBeta) {
        beta = newBeta;
    }

    public double getBeta() {
        return beta;
    }
    
    public void setN(int newN) {
    	n = newN;
    }
	
    public int getN() {
    	return n;
    }
    
    public static double dudphi(int i, double cosp, double sinp){
		
		//double cos2p = 2*cosp*cosp - 1;
		double sin2p = 2*sinp*cosp;
		
		double cos3p = 4*cosp*cosp*cosp - 3*cosp;
		double sin3p = 3*sinp - 4*sinp*sinp*sinp;
		
		//double cos6p = 2*cos3p*cos3p - 1;
		double sin6p = 2*cos3p*sin3p;
    	
		switch(i) {
		case 2: return  2 *6290.2075 *sin2p;
		case 3: return -3 *503.2166  *sin3p;
		case 6: return  6 *251.6083  *sin6p;
		
		default: throw new IllegalArgumentException("Happens only when n equals 2, 3, 6");
		}
		
    }
	
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }
	
    protected final Vector3D dr12, dr23, dr34, dra, drb, drbxa;
    protected final Vector[] gradient;
    protected double phiEq, beta;
    protected int n;
    protected Boundary boundary;
    private static final long serialVersionUID = 1L;
}
