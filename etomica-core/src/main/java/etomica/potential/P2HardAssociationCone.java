/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;
import etomica.atom.IAtomList;
import etomica.api.IBoundary;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.atom.IAtomOriented;
import etomica.space.Space;
import etomica.units.Angle;
import etomica.units.Dimension;
import etomica.units.Energy;
import etomica.units.Length;
import etomica.units.Null;

/**
 * Lennard-Jones potential with a square-well cone of attraction. 
 *
 * @author Jayant K. Singh
 */

public class P2HardAssociationCone extends Potential2 {
    private static final long serialVersionUID = 1L;
    public static boolean FLAG = false;
    private double wellcutoffFactor;
    private double wellCutoffSquared;
    private double sigma, sigmaSquared;
    private double epsilon, epsilon4, wellEpsilon;
    private double cutoffLJSquared, cutoffFactor;
    private double ec2;
    private final Vector dr;
    private IBoundary boundary;
    
    public P2HardAssociationCone(Space space) {
        this(space, 1.0, 1.0, 2.0, 8.0);
    }
    
    public P2HardAssociationCone(Space space, double sigma, double epsilon, double cutoffFactorLJ, double wellConstant) {
        super(space);
        dr = space.makeVector();

        setSigma(sigma);
        setEpsilon(epsilon);
        setCutoffFactorLJ(cutoffFactorLJ);
        setWellCutoffFactor(1.0);
        setWellEpsilon(wellConstant*getEpsilon());
        setTheta(etomica.units.Degree.UNIT.toSim(27.0));
    }

    /**
     * Returns infinity.
     */
    public double getRange() {
        return sigma*cutoffFactor;
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
                 
       // FLAG = false;
        if(r2 > cutoffLJSquared) {
            eTot = 0.0;
        }
        else {
            double s2 = sigmaSquared/r2;
            double s6 = s2*s2*s2;
            eTot = epsilon4*s6*(s6 - 1.0);
        }
                  
        if (r2 < wellCutoffSquared) {
            Vector e1 = atom0.getOrientation().getDirection();
            double er1 = e1.dot(dr);

            if ( er1 > 0.0 && er1*er1 > ec2*r2) {
                Vector e2 = atom1.getOrientation().getDirection();
                double er2 = e2.dot(dr);
                if(er2 < 0.0 && er2*er2 > ec2*r2) eTot -= wellEpsilon;
                //if(er2 < 0.0 && er2*er2 > ec2*r2) {
                	//System.out.println ("haha " + eTot);
                	//if (eTot < -2) {
                		//FLAG = true;
                	//}
                //}
                
            }
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
        setCutoffFactorLJ(cutoffFactor);
    }
    public static final Dimension getSigmaDimension() {return Length.DIMENSION;}

    /**
    * Accessor method for Lennard-Jones cutoff distance; divided by sigma
    * @return cutoff distance, divided by size parameter (sigma)
    */
    public double getCutoffFactorLJ() {return cutoffFactor;}
    /**
     * Accessor method for Lennard-Jones cutoff distance; divided by sigma
     * @param rc cutoff distance, divided by size parameter (sigma)
     */
    public void setCutoffFactorLJ(double rc) {  
        cutoffFactor = rc;
        double cutoffLJ = sigma*cutoffFactor;
        cutoffLJSquared = cutoffLJ*cutoffLJ;
    }
    public static final Dimension getCutoffFactorLJDimension() {return Null.DIMENSION;}
   
    /**
    * Accessor method for attractive-well diameter divided by sigma
    */
    public double getWellCutoffFactor() {return wellcutoffFactor;}
    /**
    * Accessor method for attractive-well diameter divided by sigma;
    */
    public void setWellCutoffFactor(double wcut) {
        wellcutoffFactor = wcut;
        double wellCutoff = sigma*wcut;
        wellCutoffSquared = wellCutoff*wellCutoff;
    }
          
    public static final Dimension getWellCutoffFactorDimension() {return Null.DIMENSION;}

    /**
    * Accessor method for Lennard-Jones energy parameter
    */ 
    public double getEpsilon() {return epsilon;}
    /**
    * Accessor method for depth of well
    */
    public void setEpsilon(double eps) {
        epsilon = eps;
        epsilon4 = 4.0 * eps;
    }
    public static final Dimension getEpsilonDimension() {return Energy.DIMENSION;}
    
    /**
    * Accessor method for attractive-well depth parameter.
    */
    public double getWellEpsilon() {return wellEpsilon;}
    /**
    * Accessor method for attractive-well depth parameter.
    */
    public void setWellEpsilon(double weps) {wellEpsilon = weps;}
          
    public static final Dimension getWellEpsilonDimension() {return Energy.DIMENSION;}
    
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
