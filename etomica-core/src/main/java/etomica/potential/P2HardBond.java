/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;
import etomica.util.Debug;

/**
 * Potential that acts like a hard string connecting the centers of two atoms.
 * Meant for use as an intra-molecular interaction. Atoms can fluctuate between
 * a minimum and maximum separation. Atoms undergo an attractive collision when
 * attempting to separate by more than the maximum and an repulsive collision
 * when attempting to come together closer than the min distance. P2Tether is 
 * similar to this potential, but does not have the repulsive core.
 */
public class P2HardBond extends Potential2HardSpherical {

    public P2HardBond(Space space) {
        this(space, 1.0, 0.15, false);
    }

    public P2HardBond(Space space, double bondLength, double bondDelta, boolean ignoreOverlap) {
        super(space);
        setBondLength(bondLength);
        setBondDelta(bondDelta);
        lastCollisionVirialTensor = space.makeTensor();
        dv = space.makeVector();
        this.ignoreOverlap = ignoreOverlap;
    }

    /**
     * Accessor method for the bond length
     */
    public final double getBondLength() {
        return bondLength;
    }

    /**
     * Accessor method for the bond extension factor
     */
    public final double getBondDelta() {
        return bondDelta;
    }

    /**
     * Setter method for the bond length
     */
    public final void setBondLength(double l) {
        bondLength = l;
        double minBondLength = bondLength - bondDelta;
        maxBondLength = bondLength + bondDelta;
        minBondLengthSquared = minBondLength * minBondLength;
        maxBondLengthSquared = maxBondLength * maxBondLength;
    }

    /**
     * Sets the bond extension factor (max = length * (1+delta))
     */
    public final void setBondDelta(double d) {
        bondDelta = d;
        setBondLength(bondLength);
    }

    public final Dimension getBondLengthDimension() {
        return Length.DIMENSION;
    }

    /**
     * Implements collision dynamics for pair attempting to separate beyond
     * tether distance
     */
    public final void bump(IAtomList pair, double falseTime) {
        IAtomKinetic atom0 = (IAtomKinetic)pair.getAtom(0);
        IAtomKinetic atom1 = (IAtomKinetic)pair.getAtom(1);
        dv.Ev1Mv2(atom1.getVelocity(), atom0.getVelocity());
        
        dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
        dr.PEa1Tv1(falseTime,dv);
        boundary.nearestImage(dr);

        double r2 = dr.squared();
        double bij = dr.dot(dv);

        if (Debug.ON && !ignoreOverlap) {
            if (bij<0.0 && Math.abs(r2 - minBondLengthSquared)/minBondLengthSquared > 1.e-9) {
                throw new RuntimeException("atoms "+pair+" not at the right distance "+r2+" "+minBondLengthSquared);
            }
            else if (bij>0.0 && Math.abs(r2 - maxBondLengthSquared)/maxBondLengthSquared > 1.e-9) {
                throw new RuntimeException("atoms "+pair+" not at the right distance "+r2+" "+maxBondLengthSquared);
            }
        }
        
        double rm0 = atom0.getType().rm();
        double rm1 = atom1.getType().rm();
        
        lastCollisionVirial = 2.0 / (rm0+rm1) * bij;
        lastCollisionVirialr2 = lastCollisionVirial / r2;
        dv.Ea1Tv1(lastCollisionVirialr2,dr);
        atom0.getVelocity().PEa1Tv1(rm0,dv);
        atom1.getVelocity().PEa1Tv1(-rm1,dv);
        atom0.getPosition().PEa1Tv1(-falseTime*rm0,dv);
        atom1.getPosition().PEa1Tv1( falseTime*rm1,dv);
    }

    public final double lastCollisionVirial() {
        return lastCollisionVirial;
    }

    public final Tensor lastCollisionVirialTensor() {
        lastCollisionVirialTensor.Ev1v2(dr, dr);
        lastCollisionVirialTensor.TE(lastCollisionVirialr2);
        return lastCollisionVirialTensor;
    }

    /**
     * Time at which two atoms will reach the end of their tether, assuming
     * free-flight kinematics
     */
    public final double collisionTime(IAtomList pair, double falseTime) {
        IAtomKinetic atom0 = (IAtomKinetic)pair.getAtom(0);
        IAtomKinetic atom1 = (IAtomKinetic)pair.getAtom(1);
        dv.Ev1Mv2(atom1.getVelocity(), atom0.getVelocity());
        
        dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
        dr.PEa1Tv1(falseTime,dv);
        boundary.nearestImage(dr);

        double r2 = dr.squared();
        double bij = dr.dot(dv);
        double v2 = dv.squared();
        
        if (ignoreOverlap && ((r2 > maxBondLengthSquared && bij > 0.0) ||
                (r2 < minBondLengthSquared && bij < 0.0))) {
            //outside bond, moving apart or overalpped and moving together; collide now
            return falseTime;
        }
        if (Debug.ON && Debug.DEBUG_NOW && ((r2 > maxBondLengthSquared && bij > 0.0) ||
                (r2 < minBondLengthSquared && bij < 0.0))) {
            System.out.println("in P2HardBond.collisionTime, "+pair+" "+r2+" "+bij+" "+maxBondLengthSquared);
            System.out.println(atom0.getPosition());
            System.out.println(atom1.getPosition());
            throw new RuntimeException("overlap");
        }
        double discr;
        if (bij < 0.0) {
            discr = bij*bij - v2 * (r2 - minBondLengthSquared);
            if(discr > 0) {
                return (-bij - Math.sqrt(discr))/v2 +falseTime;
            }
        }

        discr = bij * bij - v2 * (r2 - maxBondLengthSquared);
        if (Debug.ON && Debug.DEBUG_NOW && ((r2 > maxBondLengthSquared && bij > 0.0) ||
                (r2 < minBondLengthSquared && bij < 0.0))) {
            System.out.println("in P2HardBond.collisionTime, "+v2+" "+discr);
        }
        return (-bij + Math.sqrt(discr)) / v2 + falseTime;
    }

    /**
     * Returns 0 if the bond is within the required distance, infinity if not.
     */
    public double u(double r2) {
        if (r2 > minBondLengthSquared && r2 < maxBondLengthSquared) return 0.0;
        return Double.POSITIVE_INFINITY;
    }
    
    /**
     * Returns infinity.  
     */
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }
    
    public double energyChange() {return 0.0;}

    private static final long serialVersionUID = 1L;
    private double minBondLengthSquared;
    private double maxBondLength, maxBondLengthSquared;
    private double bondLength;
    private double bondDelta;
    private double lastCollisionVirial = 0.0;
    private double lastCollisionVirialr2 = 0.0;
    private boolean ignoreOverlap;
    private final Vector dv;
    private final Tensor lastCollisionVirialTensor;
}
