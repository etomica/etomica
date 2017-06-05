/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.atom.IAtom;
import etomica.space.Boundary;
import etomica.box.Box;
import etomica.util.random.IRandom;
import etomica.space.Vector;
import etomica.atom.IAtomOriented;
import etomica.space.Space;

public class BiasVolumeSphereOrientedDoubleSites extends BiasVolume {
    
    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private double radius;
    private double innerRadius;
    private final Vector work;
    private final Vector direction;
    private final IRandom random;
    private Boundary boundary;
    private double ec2;
    private double ec1;
    
    public BiasVolumeSphereOrientedDoubleSites(Space space, IRandom random){
        super(space);
        this.random = random;
        work = space.makeVector();
        direction = space.makeVector();
        radius = 1.0;//size of the sphere
        innerRadius = 0.8;
    }
    
    public void setBox(Box box) {
    	boundary = box.getBoundary();
    }
    
    public void setBiasSphereRadius(double radius) {
        this.radius = radius;
    }
    
    public void setBiasSphereInnerRadius(double innerRadius) {
        this.innerRadius = innerRadius;
    }
    
    public double getBiasSphereRadius() {return radius;}
    
    public double biasVolume() {
    	double partialVolume = Math.PI*(radius*radius*radius-innerRadius*innerRadius*innerRadius)*(1.0/3.0*(ec1-ec1*ec2)+2.0/3.0-ec1+1.0/3.0*ec1*ec2);//twice??
        double totalVolume = 4.0/3.0*Math.PI*(radius*radius*radius-innerRadius*innerRadius*innerRadius);
        return 2*partialVolume*partialVolume/totalVolume;//2 cones
    }
    
    
    // Insert atom1 in to the Bonding region of atom2 , randomly rotate/move atom1 and check the atom2 and then check the atom1
    //Bonding region is a cube /rectangle for 2 D
    public void biasInsert(IAtom atom1, IAtom atom2) {
        double er2;
        do {
            work.setRandomInSphere(random);//work =bond vector from atom2 to atom1
            work.TE(radius);//multiply by radius
			while (work.squared()<innerRadius*innerRadius){
				work.setRandomInSphere(random);
				work.TE(radius);
			}
			//compute the orientation
			Vector e2 = ((IAtomOriented)atom2).getOrientation().getDirection();//orientation of atom2
			er2 = e2.dot(work);
        }
        while ( er2*er2 < ec2*work.squared());
        atom1.getPosition().Ev1Pv2(atom2.getPosition(), work);//changing the position of atom1(atom2 is fixed)
        //rotate atom2 to have a right orientation
        double er1;
        do{
        	direction.setRandomSphere(random);//finding the orientation of atom1 randomly
            er1 = direction.dot(work);
        }
        while ( er1*er2 < 0.0 || er1*er1 < ec2*work.squared());//two atoms are bonding
        ((IAtomOriented)atom1).getOrientation().setDirection(direction);//once we find the one, then set the orientation of atom1(atom2 is fixed)
        
    }

    /**
     *Function to check for bonding
     *calculate the distance between atom1 and atom2
     */
    
    public boolean isAssociated(IAtom atom1, IAtom atom2){
    
        work.E(atom2.getPosition());
        work.ME(atom1.getPosition());
        boundary.nearestImage(work);
        double r2 = work.squared();
//        if (atom1.getLeafIndex() == 68 ||atom2.getLeafIndex() == 68 ||atom1.getLeafIndex() == 303 ||atom2.getLeafIndex() == 303){
//	        System.out.println ("atom1 = "+atom1+ " atom2 = "+atom2);
//	        System.out.println ("r2 = "+r2);
//        }
        if (r2 < innerRadius*innerRadius || r2 > radius*radius) {
        	return false;
        }
        Vector e1 = ((IAtomOriented)atom1).getOrientation().getDirection();
        double er1 = e1.dot(work);
//        if (atom1.getLeafIndex() == 68 ||atom2.getLeafIndex() == 68 ||atom1.getLeafIndex() == 303 ||atom2.getLeafIndex() == 303){
//	        System.out.println ("atom1 = "+atom1+ " atom2 = "+atom2);
//	        System.out.println ("er1 = "+er1);
//        }
        if ( er1*er1 < ec2*r2) {
        	return false;
        }
        Vector e2 = ((IAtomOriented)atom2).getOrientation().getDirection();
        double er2 = e2.dot(work);
//        if (atom1.getLeafIndex() == 68 ||atom2.getLeafIndex() == 68 ||atom1.getLeafIndex() == 303 ||atom2.getLeafIndex() == 303){
//	        System.out.println ("atom1 = "+atom1+ " atom2 = "+atom2);
//	        System.out.println ("er2 = "+er2);
//        }       
        return er1*er2 > 0.0 && er2*er2 > ec2*r2;   
    }
public double getTheta() {return Math.acos(ec1);}
    
    /**
     * Accessor method for angle (in radians) describing width of cone.
     */
    public void setTheta(double t) {
        ec1   = Math.cos(t);
        ec2   = ec1*ec1;
    }
}
