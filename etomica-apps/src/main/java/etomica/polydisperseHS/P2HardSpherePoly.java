/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.polydisperseHS;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.Potential2;
import etomica.potential.Potential2HardSpherical;
import etomica.potential.PotentialHard;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;
import etomica.util.Debug;

import java.util.Arrays;

/**
 * Basic hard-(rod/disk/sphere) potential.
 * Energy is infinite if spheres overlap, and is zero otherwise.  Collision diameter describes
 * size of spheres.
 * Suitable for use in space of any dimension.
 *
 * @author David Kofke
 */
public class P2HardSpherePoly extends Potential2 implements PotentialHard {

    protected double lastCollisionVirial = 0.0;
    protected double lastCollisionVirialr2 = 0.0;
    protected final boolean ignoreOverlap;
    protected final Vector dr, dv;
    protected final double[] sigmas;
    protected final Tensor lastCollisionVirialTensor;
    protected Boundary boundary;

    public P2HardSpherePoly(Space space, double[] sigmas, boolean ignoreOverlap) {
        super(space);
        this.sigmas = sigmas;
        lastCollisionVirialTensor = space.makeTensor();
        dr = space.makeVector();
        dv = space.makeVector();
        this.ignoreOverlap = ignoreOverlap;
    }

    public double getRange() {
        double max = sigmas[0];
        for(int i=1;i<sigmas.length;i++){
            if(sigmas[i]>max){
                max = sigmas[i];
            }
        }
        return max;
    }

    /**
     * Time to collision of pair, assuming free-flight kinematics
     */
    public double collisionTime(IAtomList pair, double falseTime) {
        IAtomKinetic atom0 = (IAtomKinetic)pair.get(0);
        IAtomKinetic atom1 = (IAtomKinetic)pair.get(1);
        dv.Ev1Mv2(atom1.getVelocity(), atom0.getVelocity());


        dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
        dr.PEa1Tv1(falseTime,dv);
        boundary.nearestImage(dr);

        double bij = dr.dot(dv);
        double time = Double.POSITIVE_INFINITY;

        if(bij < 0.0) {
            double v2 = dv.squared();
            double sig2 = (sigmas[atom0.getLeafIndex()]+sigmas[atom1.getLeafIndex()])/2;
            sig2 *= sig2;
            if (ignoreOverlap && dr.squared() < sig2) return falseTime+0.001*Math.sqrt(dr.squared())/Math.sqrt(v2);
            double discriminant = bij*bij - v2 * ( dr.squared() - sig2 );
            if(discriminant > 0) {
                time = (-bij - Math.sqrt(discriminant))/v2;
            }
        }
        if (Debug.ON && Debug.DEBUG_NOW && (Debug.allAtoms(pair) || time < 0.0)) {
        	System.out.println("atoms "+pair+" r2 "+dr.squared()+" bij "+bij+" time "+time);
        	if (time < 0.0) throw new RuntimeException("negative collision time for hard spheres");
        }
        return time + falseTime;
    }
    
    /**
     * Implements collision dynamics and updates lastCollisionVirial
     */
    public void bump(IAtomList pair, double falseTime) {
        IAtomKinetic atom0 = (IAtomKinetic)pair.get(0);
        IAtomKinetic atom1 = (IAtomKinetic)pair.get(1);
        dv.Ev1Mv2(atom1.getVelocity(), atom0.getVelocity());
        
        dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
        dr.PEa1Tv1(falseTime,dv);
        boundary.nearestImage(dr);

        double r2 = dr.squared();
        double bij = dr.dot(dv);
        double rm0 = atom0.getType().rm();
        double rm1 = atom1.getType().rm();
        double reducedMass = 2.0/(rm0 + rm1);
        lastCollisionVirial = reducedMass*bij;
        lastCollisionVirialr2 = lastCollisionVirial/r2;
        //dv is now change in velocity due to collision
        dv.Ea1Tv1(lastCollisionVirialr2,dr);
        atom0.getVelocity().PEa1Tv1( rm0,dv);
        atom1.getVelocity().PEa1Tv1(-rm1,dv);
        atom0.getPosition().PEa1Tv1(-falseTime*rm0,dv);
        atom1.getPosition().PEa1Tv1( falseTime*rm1,dv);
    }
    
    public double lastCollisionVirial() {
        return lastCollisionVirial;
    }
    
    public Tensor lastCollisionVirialTensor() {
        lastCollisionVirialTensor.Ev1v2(dr, dr);
//        lastCollisionVirialTensor.DE(lastCollisionVirialr2);
        return lastCollisionVirialTensor;        
    }
    
    /**
     * Accessor method for collision diameter
     */
    public double[] getCollisionDiameter() {return sigmas;}

    /**
     * Interaction energy of the pair.
     * Zero if separation is greater than collision diameter, infinity otherwise
     */
    public double energy(IAtomList pair) {
        IAtom atom0 = pair.get(0);
        IAtom atom1 = pair.get(1);

        dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
        boundary.nearestImage(dr);

        double sig2 = (sigmas[atom0.getLeafIndex()]+sigmas[atom1.getLeafIndex()])/2;
        sig2 *= sig2;

        if (dr.squared() < sig2) {
            return Double.POSITIVE_INFINITY;
        }
        return 0.0;

    }

    public void setBox(Box box) {
        boundary = box.getBoundary();
    }




    public double energyChange() {return 0.0;}
    
}//end of P2HardSphere
