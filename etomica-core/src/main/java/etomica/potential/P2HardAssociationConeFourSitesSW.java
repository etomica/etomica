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
import etomica.space3d.IOrientationFull3D;
import etomica.units.dimensions.Angle;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;

/**
 * Square-Well cone of attraction for associating fluid theory
 * Four association sites
 *
 * @author Hye Min Kim
 */

public class P2HardAssociationConeFourSitesSW extends Potential2 {
    private static final long serialVersionUID = 1L;
    public static boolean FLAG = false;
    private double wellCutoffFactor, innerWellCutoffFactor;
    private double wellCutoffSquared, innerWellCutoffSquared;
    private double sigma, sigmaSquared;
    private double epsilon, epsilon4, wellEpsilon;
    private double cutoffLJSquared, cutoffFactor;
    private double ec1, ec2;
    private final Vector dr;
    private Boundary boundary;
    
    public P2HardAssociationConeFourSitesSW(Space space, double sigma, double epsilon, double cutoffFactorLJ, double wellConstant) {
        super(space);
        dr = space.makeVector();

        setSigma(sigma);
        setEpsilon(epsilon);
        setCutoffFactorLJ(cutoffFactorLJ);//truncation distance of potential
        setWellCutoffFactor(1.0);
        setInnerWellCutoffFactor(0.0);
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
        IAtomOriented atom0 = (IAtomOriented)atoms.get(0);
        IAtomOriented atom1 = (IAtomOriented)atoms.get(1);
        dr.Ev1Mv2(atom1.getPosition(),atom0.getPosition());//dr=atom1-atom0
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        double eTot = 0.0;
                 
       // FLAG = false;
        if(r2 > wellCutoffSquared || r2 < innerWellCutoffSquared) {
            eTot = 0.0;
        }                 
        else {
        	double angle1 = Math.acos(-1.0/3.0);//theta = 109.471Math.acos(-1.0/3.0)
        	double cosAngle1 = -1.0/3.0;
        	double sinAngle1 = Math.sin(angle1);
        	//double angle2 = (Math.PI-angle1)/2.0;//(pi-theta)/2
        	//double cosAngle2 = Math.cos(angle2);
        	//double length = 2.0*cosAngle2;//length of side of tetrahedron
        	double coordX = cosAngle1;//X=-0.333
            double coordY = cosAngle1*(1-cosAngle1)/sinAngle1;//Y=-0.471
            double coordZ = Math.sqrt(1-coordX*coordX-coordY*coordY);//Z=0.816
//            System.out.println("coorX: "+coordX);
//            System.out.println("coorY: "+coordY);
//            System.out.println("coorZ: "+coordZ);
            //double coordZ = 0.5*length;
        	Vector e1A = atom0.getOrientation().getDirection();
            double er1Aa = e1A.dot(dr);//vector of site Aa on atom0  
            Vector e1AaY = ((IOrientationFull3D) atom0.getOrientation()).getSecondaryDirection();//Perpendicular(second) direction of e1A
            Vector e1AaZ = space.makeVector();//third direction of e1A
            e1AaZ.E(e1A);
            e1AaZ.XE(e1AaY);//crossproduct of e1A and e1AaY
            Vector e1Ab = space.makeVector();
            e1Ab.Ea1Tv1(cosAngle1, e1A);
            e1Ab.PEa1Tv1(sinAngle1, e1AaY);
            double er1Ab = e1Ab.dot(dr);//vector of site Ab on atom0
            Vector e1Ba = space.makeVector();
            e1Ba.Ea1Tv1(coordX, e1A);
            e1Ba.PEa1Tv1(coordY, e1AaY);
            e1Ba.PEa1Tv1(coordZ, e1AaZ);
            double er1Ba = e1Ba.dot(dr);//vector of site Ba on atom0
            Vector e1Bb = space.makeVector();
            e1Bb.Ea1Tv1(coordX, e1A);
            e1Bb.PEa1Tv1(coordY, e1AaY);
            e1Bb.PEa1Tv1(-coordZ, e1AaZ);
            double er1Bb = e1Bb.dot(dr);//vector of site Bb on atom0
            Vector e2A = atom1.getOrientation().getDirection();
            double er2Aa = e2A.dot(dr);//vector of site Aa on atom1
            Vector e2AaY = ((IOrientationFull3D) atom1.getOrientation()).getSecondaryDirection();//Perpendicular direction of e2A
            Vector e2AaZ = space.makeVector();
            e2AaZ.E(e2A);
            e2AaZ.XE(e2AaY);//crossproduct
            Vector e2Ab = space.makeVector();
            e2Ab.Ea1Tv1(cosAngle1, e2A);
            e2Ab.PEa1Tv1(sinAngle1, e2AaY);
            double er2Ab = e2Ab.dot(dr);//vector of site Ab on atom1
            Vector e2Ba = space.makeVector();
            e2Ba.Ea1Tv1(coordX, e2A);
            e2Ba.PEa1Tv1(coordY, e2AaY);
            e2Ba.PEa1Tv1(coordZ, e2AaZ);
            double er2Ba = e2Ba.dot(dr);//vector of site Ba on atom1
            Vector e2Bb = space.makeVector();
            e2Bb.Ea1Tv1(coordX, e2A);
            e2Bb.PEa1Tv1(coordY, e2AaY);
            e2Bb.PEa1Tv1(-coordZ, e2AaZ);
            double er2Bb = e2Bb.dot(dr);//vector of site Bb on atom1
            
            if (er1Aa*er2Aa < 0.0 ||er1Aa*er2Ab < 0.0||er1Ab*er2Aa < 0.0||er1Ab*er2Ab < 0.0||er1Ba*er2Ba < 0.0 ||er1Ba*er2Bb < 0.0||er1Bb*er2Ba < 0.0||er1Bb*er2Bb < 0.0 ){
            	eTot -= 0.0; 
            }
            
            if ((er1Aa*er2Ba < 0.0 && er1Aa*er1Aa > ec2*r2 && er2Ba*er2Ba > ec2*r2) || (er1Aa*er2Bb < 0.0 && er1Aa*er1Aa > ec2*r2 && er2Bb*er2Bb > ec2*r2)
                ||(er1Ab*er2Ba < 0.0 && er1Ab*er1Ab > ec2*r2 && er2Ba*er2Ba > ec2*r2)||(er1Ab*er2Bb < 0.0 && er1Ab*er1Ab > ec2*r2 && er2Bb*er2Bb > ec2*r2)
                ||(er1Ba*er2Aa < 0.0 && er1Ba*er1Ba > ec2*r2 && er2Aa*er2Aa > ec2*r2)|| (er1Ba*er2Ab < 0.0 && er1Ba*er1Ba > ec2*r2 && er2Ab*er2Ab > ec2*r2)
                ||(er1Bb*er2Aa < 0.0 && er1Bb*er1Bb > ec2*r2 && er2Aa*er2Aa > ec2*r2)||(er1Bb*er2Ab < 0.0 && er1Bb*er1Bb > ec2*r2 && er2Ab*er2Ab > ec2*r2)) 
            {
            	eTot -= wellEpsilon;
            	//System.out.println("eTot= "+eTot);
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
    public double getWellCutoffFactor() {return wellCutoffFactor;}
    /**
    * Accessor method for attractive-well diameter divided by sigma;
    */
    public void setWellCutoffFactor(double wcut) {
        wellCutoffFactor = wcut;
        double wellCutoff = sigma*wcut;
        wellCutoffSquared = wellCutoff*wellCutoff;
    }
    
    public void setInnerWellCutoffFactor(double wcut) {
        innerWellCutoffFactor = wcut;
        double innerWellCutoff = sigma*wcut;
        innerWellCutoffSquared = innerWellCutoff*innerWellCutoff;
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
    public double getTheta() {return Math.acos(ec1);}
    
    /**
     * Accessor method for angle (in radians) describing width of cone.
     */
    public void setTheta(double t) {
        ec1   = Math.cos(t);
        ec2   = ec1*ec1;
    }
    public Dimension getThetaDimension() {return Angle.DIMENSION;}

    public void setBox(Box box) {
        boundary = box.getBoundary();
    }
}
