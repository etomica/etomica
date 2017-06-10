/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;
import etomica.atom.IAtomList;
import etomica.space.Boundary;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.atom.IAtomOriented;
import etomica.space.Space;
import etomica.units.dimensions.Angle;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;

/**
 * potential with a square-well cone of attraction for reference system. 
 *
 *
 * @author Hye Min Kim
 */

public class P2HardAssociationConeReference extends Potential2 {
    private static final long serialVersionUID = 1L;
    public static boolean FLAG = false;
    private double sigma, sigmaSquared;
    private double ec2;
    private final Vector dr;
    private Boundary boundary;
    
    public P2HardAssociationConeReference(Space space, double sigma) {
        super(space);
        dr = space.makeVector();

        setSigma(sigma);
        setTheta(etomica.units.Degree.UNIT.toSim(27.0));
    }

    /**
     * Returns infinity.
     */
    public double getRange() {
        return sigma;
    }


    /**
     * Returns the pair potential energy.
     */
    public double energy(IAtomList atoms) {
        IAtomOriented atom0 = (IAtomOriented)atoms.getAtom(0);
        IAtomOriented atom1 = (IAtomOriented)atoms.getAtom(1);
        dr.Ev1Mv2(atom1.getPosition(),atom0.getPosition());
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        double eTot = 0.0;
                
                  
        if (r2 < sigmaSquared) {
        	Vector e1A = atom0.getOrientation().getDirection();
            double er1A = e1A.dot(dr);//vector of site A on atom0
            double er1B = -1*er1A;//vector of site B on atom0
            Vector e2A = atom1.getOrientation().getDirection();
            double er2A = e2A.dot(dr);//vector of site A on atom1
            double er2B = -1*er2A;//vector of site B on atom1
            if(er1A*er2B < 0.0 && er1A*er1A > ec2*r2 && er2B*er2B > ec2*r2) eTot = Double.POSITIVE_INFINITY;
                //if(er2 < 0.0 && er2*er2 > ec2*r2) {
                	//System.out.println ("haha " + eTot);
                	//if (eTot < -2) {
                		//FLAG = true;
                	//}
                //}
        } 
        return eTot;
    }
    
    /**
     * Accessor method for Lennard-Jones size parameter
     */
    public double getSigma() {return sigma;}
    /**
     * Accessor method for Lennard-Jones size parameter
     */
    public void setSigma(double s) {
        sigma = s;
        sigmaSquared = s*s;
    }
    public static final Dimension getSigmaDimension() {return Length.DIMENSION;}

    /**
     * Accessor method for angle describing width of cone.
     */
    public double getTheta() {return Math.acos(ec2);}
    
    /**
     * Accessor method for angle (in radians) describing width of cone.
     */
    public void setTheta(double t) {
        ec2   = Math.cos(t);
        ec2   = ec2*ec2;
    }
    public Dimension getThetaDimension() {return Angle.DIMENSION;}

    public void setBox(Box box) {
        boundary = box.getBoundary();
    }

}
