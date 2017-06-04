/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.atom.IAtomOrientedKinetic;
import etomica.atom.IAtomTypeOriented;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space3d.Vector3D;

/**
 * Basic hard-(rod/disk/sphere) potential, with surface roughness to couple rotation and translational motions.
 * Suitable for use in space of any dimension.
 *
 * @author David Kofke
 */
public class P2RoughSphere extends P2HardSphere {

    private static final long serialVersionUID = 1L;
    private final Vector3D omegaSum = new Vector3D();
    private final Vector3D v12Surface;
    private final Vector3D v12Par;
    private final Vector3D v12Perp;
    private final Vector3D impulse;
    
    public P2RoughSphere(Space space) {
        this(space, 1.0, false);
    }
    
    public P2RoughSphere(Space space, double d, boolean ignoreOverlap) {
        super(space,d,ignoreOverlap);
        v12Surface = new Vector3D();
        v12Par = new Vector3D();
        v12Perp = new Vector3D();
        impulse = new Vector3D();
    }

    /**
     * Implements collision dynamics and updates lastCollisionVirial
     * Assumes atoms have same size and mass
     */
    public void bump(IAtomList pair, double falseTime) {
        IAtom atom0 = pair.getAtom(0);
        IAtom atom1 = pair.getAtom(1);
        IAtomOrientedKinetic coord0 = (IAtomOrientedKinetic)atom0;
        IAtomOrientedKinetic coord1 = (IAtomOrientedKinetic)atom1;
        Vector v1 = coord0.getVelocity();
        Vector v2 = coord1.getVelocity();
        dv.Ev1Mv2(v2, v1);
        
        dr.Ev1Mv2(coord1.getPosition(), coord0.getPosition());
        dr.PEa1Tv1(falseTime,dv);
        boundary.nearestImage(dr);

        double r2 = dr.squared();
        double bij = dr.dot(dv);
        double rm0 = atom0.getType().rm();
        double rm1 = atom1.getType().rm();
        double kappa = 4*((IAtomTypeOriented)atom0.getType()).getMomentOfInertia().getX(0)*rm0/(collisionDiameter*collisionDiameter);
        omegaSum.E(coord0.getAngularVelocity());
        omegaSum.PE(coord1.getAngularVelocity());
        // v12Surface should come to equal v2 - v1 - 1/2*(omega2+omega1) X (r2-r1)
        v12Surface.E(dr); // (r2 - r1)
        v12Surface.XE(omegaSum); //(r2-r1) X (omega2+omega1)
        v12Surface.TE(0.5); // +1/2 (r2-r1) X (omega2+omega1) [which equals -1/2*(omega2+omega1) X (r2-r1)]
        v12Surface.PE(v2);// p2/m2 +1/2 (r2-r1) X (omega2+omega1)
        v12Surface.ME(v1);// p2/m2 - p1/m1 +1/2 (r2-r1) X (omega2+omega1)
        //component of v12Surface parallel to r2-r1: v12Par = (v12Surface . dr) dr / |dr|^2
        v12Par.E(dr);
        v12Par.TE(v12Surface.dot(dr)/r2);
        //component of v12Surface perpendicular to r2-r1:  v12Perp = v12Surface - v12Par
        v12Perp.E(v12Surface);
        v12Perp.ME(v12Par);
        
        impulse.E(v12Par);
        impulse.PEa1Tv1(kappa/(1+kappa),v12Perp);
        impulse.TE(atom0.getType().getMass());
        
        coord0.getVelocity().PEa1Tv1( rm0,impulse);
        coord1.getVelocity().PEa1Tv1(-rm1,impulse);
        coord0.getPosition().PEa1Tv1(-falseTime*rm0,impulse);
        coord1.getPosition().PEa1Tv1( falseTime*rm1,impulse);
        
        //here omegaSum is used to hold the angular impulse
        omegaSum.E(dr);
        omegaSum.XE(impulse);
        omegaSum.TE(-0.5);
        coord0.getAngularVelocity().PE(omegaSum);
        coord1.getAngularVelocity().PE(omegaSum);
        
        lastCollisionVirial = 2.0/(rm0 + rm1)*bij;
        lastCollisionVirialr2 = lastCollisionVirial/r2;
    }
    
    //need to consider if hard-sphere virial is same as rough sphere virial
    public final double lastCollisionVirial() {
        return Double.NaN;
      //  return lastCollisionVirial;
    }
    
    //especially need to consider more carefully this method
    public final Tensor lastCollisionVirialTensor() {
        lastCollisionVirialTensor.Ev1v2(dr, dr);
        lastCollisionVirialTensor.TE(lastCollisionVirialr2);
        return lastCollisionVirialTensor;        
    }
}
