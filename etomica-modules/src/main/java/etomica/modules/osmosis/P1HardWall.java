/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.osmosis;

import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.potential.Potential1;
import etomica.potential.PotentialHard;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;

/**
 */
 
public class P1HardWall extends Potential1 implements PotentialHard {
    
    private static final long serialVersionUID = 1L;
    private double collisionRadius;
    
    public P1HardWall(Space space) {
        this(space, 1.0);
    }
    
    public P1HardWall(Space space, double sigma) {
        super(space);
        collisionRadius = sigma;
    }

    public double energy(IAtomList a) {
        double e = 0.0;
        //XXX ignore atoms in the wall.  this can happen due to bogus initial configurations
//        if (Math.abs(((AtomLeaf)a).coord.position().get(0)) < collisionRadius) {
//            e = Double.MAX_VALUE;
//        }
        return e;
    }

     
    public double collisionTime(IAtomList a, double falseTime) {
        IAtomKinetic atom = (IAtomKinetic)a.getAtom(0);
        Vector r = atom.getPosition();
        Vector v = atom.getVelocity();
        double vx = v.getX(0);
        double rx = r.getX(0) + vx * falseTime;
        double t = (vx > 0.0) ? - collisionRadius : collisionRadius;
        t = (t - rx) / vx;
        if (t < 0) {
            // moving away from the wall
            t = Double.POSITIVE_INFINITY;
        }
        return t+falseTime;
    }

    public void bump(IAtomList a, double falseTime) {
        IAtomKinetic atom = (IAtomKinetic)a.getAtom(0);
        Vector v = atom.getVelocity();

        v.setX(0,-v.getX(0));

        double newP = atom.getPosition().getX(0) - falseTime*v.getX(0)*2.0;
        atom.getPosition().setX(0,newP);
    }

    public double energyChange() {
        return 0;
    }
    
    /**
     * not yet implemented
     */
    public double lastCollisionVirial() {return Double.NaN;}
    
    /**
     * not yet implemented.
     */
    public Tensor lastCollisionVirialTensor() {return null;}
    
        
    /**
     * Distance from the center of the sphere to the boundary at collision.
     */
    public void setCollisionRadius(double d) {collisionRadius = d;}
    /**
     * Distance from the center of the sphere to the boundary at collision.
     */
    public double getCollisionRadius() {return collisionRadius;}
    /**
     * Indicates collision radius has dimensions of Length.
     */
    public Dimension getCollisionRadiusDimension() {return Length.DIMENSION;}

}
   
