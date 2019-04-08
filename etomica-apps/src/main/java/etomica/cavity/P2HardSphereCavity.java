/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.cavity;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.potential.P2HardSphere;
import etomica.space.Space;
import etomica.util.random.IRandom;

public class P2HardSphereCavity extends P2HardSphere {

    protected final IRandom random;
    protected double pairWell;
    protected IAtomKinetic pairedAtom1, pairedAtom2;
    protected int idxSum, idxProduct;

    public enum CollisionType {INTERNAL_BOUNCE, ESCAPE, EXTERNAL_BOUNCE, CAPTURE}

    protected CollisionType lastCollisionType;

    public P2HardSphereCavity(Space space) {
        this(space, null, 0);
    }

    public P2HardSphereCavity(Space space, IRandom random, double pairWell) {
        super(space, 1, false);
        this.random = random;
        this.pairWell = pairWell;
    }

    public IAtom getPairedAtom1() {
        return pairedAtom1;
    }

    public IAtom getPairedAtom2() {
        return pairedAtom2;
    }

    public double collisionTime(IAtomList pair, double falseTime) {
        int idx0 = pair.get(0).getLeafIndex();
        int idx1 = pair.get(1).getLeafIndex();
        boolean paired = idx0 * idx1 == idxProduct && idx0 + idx1 == idxSum;
        if (!paired) {
            // this pair is not overlapping, treat as normal
            return super.collisionTime(pair, falseTime);
        }
        // this pair is overlapping.  look for escape
        IAtomKinetic atom0 = (IAtomKinetic) pair.get(0);
        IAtomKinetic atom1 = (IAtomKinetic) pair.get(1);
        dv.Ev1Mv2(atom1.getVelocity(), atom0.getVelocity());

        dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
        dr.PEa1Tv1(falseTime, dv);
        boundary.nearestImage(dr);

        double bij = dr.dot(dv);
        double v2 = dv.squared();
        double discriminant = bij * bij - v2 * (dr.squared() - sig2);
        if (discriminant < 0) {
            throw new RuntimeException("paired atoms not overlapping?");
        }

        double time = (-bij + Math.sqrt(discriminant)) / v2;
        return time + falseTime;
    }

    /**
     * Implements collision dynamics and updates lastCollisionVirial
     */
    public void bump(IAtomList pair, double falseTime) {
        int idx0 = pair.get(0).getLeafIndex();
        int idx1 = pair.get(1).getLeafIndex();
        boolean paired = idx0 * idx1 == idxProduct && idx0 + idx1 == idxSum;
        if (idxSum > 0 && (!paired || pairWell == Double.POSITIVE_INFINITY)) {
            lastCollisionType = paired ? CollisionType.INTERNAL_BOUNCE : CollisionType.EXTERNAL_BOUNCE;
            super.bump(pair, falseTime);
            return;
        }

        IAtomKinetic atom0 = (IAtomKinetic) pair.get(0);
        IAtomKinetic atom1 = (IAtomKinetic) pair.get(1);
        dv.Ev1Mv2(atom1.getVelocity(), atom0.getVelocity());

        dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
        dr.PEa1Tv1(falseTime, dv);
        boundary.nearestImage(dr);

        double r2 = dr.squared();
        double bij = dr.dot(dv);

        if (paired != (bij > 0.0)) {         // Separating
            throw new RuntimeException("paired: " + paired + " and bij " + bij);
        }

        double rm0 = atom0.getType().rm();
        double rm1 = atom1.getType().rm();
        double reducedMass = 2.0 / (rm0 + rm1);

        if (bij > 0) {
            if (pairWell > 0 && random.nextDouble() < Math.exp(-pairWell)) {     // Not enough kinetic energy to escape
                lastCollisionVirial = reducedMass * bij;
                lastCollisionType = CollisionType.INTERNAL_BOUNCE;
            } else {                 // Escape
                lastCollisionVirial = 0;
                pairedAtom1 = pairedAtom2 = null;
                idxSum = idxProduct = 0;
                lastCollisionType = CollisionType.ESCAPE;
            }
        } else if (pairWell > 0 || random.nextDouble() < Math.exp(pairWell)) {   // Approach/capture
            // core collision -- capture
            lastCollisionVirial = 0;
            pairedAtom1 = atom0;
            pairedAtom2 = atom1;
            idxSum = idx0 + idx1;
            idxProduct = idx0 * idx1;
            lastCollisionType = CollisionType.CAPTURE;
        } else {
            // not enough kinetic energy to overcome shoulder
            lastCollisionVirial = reducedMass * bij;
            lastCollisionType = CollisionType.EXTERNAL_BOUNCE;
        }

        lastCollisionVirialr2 = lastCollisionVirial / r2;
        //dv is now change in velocity due to collision
        dv.Ea1Tv1(lastCollisionVirialr2, dr);
        atom0.getVelocity().PEa1Tv1(rm0, dv);
        atom1.getVelocity().PEa1Tv1(-rm1, dv);
        atom0.getPosition().PEa1Tv1(-falseTime * rm0, dv);
        atom1.getPosition().PEa1Tv1(falseTime * rm1, dv);
    }

    public CollisionType getLastCollisionType() {
        return lastCollisionType;
    }
}
