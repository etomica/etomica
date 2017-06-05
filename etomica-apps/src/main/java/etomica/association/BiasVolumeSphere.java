/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.atom.IAtom;
import etomica.space.Boundary;
import etomica.space.Vector;
import etomica.box.Box;
import etomica.util.random.IRandom;
import etomica.space.Space;

public class BiasVolumeSphere extends BiasVolume {
    
    private static final long serialVersionUID = 1L;
    private double radius;
    private double innerRadius;
    private final Vector work;
    private final IRandom random;
    private Boundary boundary;
    
    public BiasVolumeSphere(Space space, IRandom random){
        super(space);
        this.random = random;
        work = space.makeVector();
        radius = 1.0;//size of the sphere
        innerRadius = 0.9;
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
        double prod = 4.0/3.0*Math.PI*(radius*radius*radius-innerRadius*innerRadius*innerRadius);
        return prod;
    }
    
    
    // Insert atom1 in to the Bonding region of atom2 
    //Bonding region is a cube /rectangle for 2 D
    public void biasInsert(IAtom atom1, IAtom atom2) {
        work.setRandomInSphere(random);
        work.TE(radius);//multiply by radius
        while (work.squared()<innerRadius*innerRadius){
        	work.setRandomInSphere(random);
        	work.TE(radius);
        }
        atom1.getPosition().Ev1Pv2(atom2.getPosition(), work);
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
        //System.out.println ("atom1 = "+atom1+ " atom2 = "+atom2);
        //System.out.println ("r2 = "+r2);
        return r2 > innerRadius*innerRadius && r2 < radius*radius;
    }
}
