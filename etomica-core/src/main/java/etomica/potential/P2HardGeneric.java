/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.space.Vector;

/**
 * Potential class for any hard spherical potential.  The steps in the potential
 * are defined at construction.  This implements Potential2Soft so that it can be
 * used in a PotentialCompute to get the energy.
 */
public class P2HardGeneric implements IPotential2 {

    /**
     * Setting fixOverlapDefault to true will cause any P2HardGeneric to be
     * created with fixOverlap on by default.  With fixOverlap enabled, the
     * energy of overlapping pairs will be returned as the energy of pairs
     * that are barely not overlapping and collision times will be returned
     * to help push the particles apart.
     */
    public static boolean fixOverlapDefault;

    /**
     * The pair distances where collisions should happen.
     * Distances are listed from smallest to largest.
     */
    protected final double[] collisionDistances2;

    /**
     * The pair energies where a pair of atoms is between the various
     * distances returned by getCollisionDistances().
     * <p>
     * energy[0] = u when r < collisionDistance[0]
     * energy[1] = u when collisionDistance[0] < r < collisionDistance[1]
     * energy[2] = u when collisionDistance[1] < r < collisionDistance[2]
     * ...
     * energy[n] = u when r > collisionDistance[n-1]
     * <p>
     * The last energy should always be 0.
     */
    protected final double[] energies;

    // only used to support the energy(IAtomList) legacy method

    protected final boolean fixOverlap;

    public P2HardGeneric(double[] collisionDistances, double[] energies) {
        this(collisionDistances, energies, fixOverlapDefault);
    }

    public P2HardGeneric(double[] collisionDistances, double[] energies, boolean fixOverlap) {
        if (collisionDistances.length != energies.length) {
            throw new IllegalArgumentException("collision distances and energies must have the same length");
        }
        this.collisionDistances2 = new double[collisionDistances.length];
        for (int i = 0; i < collisionDistances.length; i++) {
            collisionDistances2[i] = collisionDistances[i] * collisionDistances[i];
        }
        this.energies = energies;
        this.fixOverlap = fixOverlap;
    }

    public double getCollisionDiameter(int i) {
        return Math.sqrt(collisionDistances2[i]);
    }

    public void setCollisionDiameter(int i, double collisionDiameter) {
        collisionDistances2[i] = collisionDiameter * collisionDiameter;
    }

    public void setEnergyForState(int i, double u) {
        energies[i] = u;
    }

    /**
     * Returns the state of the pair of atoms (atom1 and atom2) at distance r12
     */
    public int getState(IAtom atom1, IAtom atom2, Vector r12) {
        double r2 = r12.squared();
        double[] cd2 = getCollisionDistances2(atom1, atom2);
        for (int i = 0; i < cd2.length; i++) {
            if (cd2[i] > r2) {
                if (i == 0 && fixOverlap) return cd2.length > 1 ? 1 : -1;
                return i;
            }
        }
        return -1;
    }

    // subclasses can override this method to have different distances for special atom pairs
    protected double[] getCollisionDistances2(IAtom atom1, IAtom atom2) {
        return collisionDistances2;
    }

    public double collisionTime(IAtomKinetic atom1, IAtomKinetic atom2, Vector r12, Vector v12, int collisionState, double falseTime) {
        return collisionTime(r12, v12, collisionState, getCollisionDistances2(atom1, atom2));
    }

    protected double collisionTime(Vector r12, Vector v12, int collisionState, double[] cd2) {

        double bij = r12.dot(v12);
        double time = Double.POSITIVE_INFINITY;

        if (collisionState < 0) collisionState = cd2.length;
        if (bij < 0.0) {
            // moving together
            double r2 = r12.squared();
            double v2 = v12.squared();
            if (fixOverlap && r2 < cd2[0]) {
                // overlapped collide now
                return 0.001 * Math.sqrt(r2 / v2);
            }
            double discriminant = collisionState == 0 ? -1 : bij * bij - v2 * (r2 - cd2[collisionState - 1]);
            if (discriminant > 0) {
                // hit
                time = (-bij - Math.sqrt(discriminant)) / v2;
            } else if (collisionState < cd2.length) {
                // miss, look for escape
                double discr = bij * bij - v2 * (r2 - cd2[collisionState]);
                time = (-bij + Math.sqrt(discr)) / v2;
            }
        } else if (collisionState < cd2.length) {
            // moving away, look for escape
            double r2 = r12.squared();
            double v2 = v12.squared();
            double discr = bij * bij - v2 * (r2 - cd2[collisionState]);
            time = (-bij + Math.sqrt(discr)) / v2;
        }
        return time;
    }

    public int bump(IAtomKinetic atom1, IAtomKinetic atom2, int oldState, Vector r12, Vector v12, double falseTime, double[] virial, double[] du) {
        double r2 = r12.squared();
        double bij = r12.dot(v12);
        double rm0 = atom1.getType().rm();
        double rm1 = atom2.getType().rm();
        double reducedMass = 1.0 / (rm0 + rm1);
        double[] cd2 = getCollisionDistances2(atom1, atom2);
        boolean core;
        double rCollision2;
        if (oldState < 0) oldState = cd2.length;
        if (oldState == cd2.length) {
            core = true;
            rCollision2 = cd2[cd2.length - 1];
        } else if (oldState == 0) {
            core = false;
            rCollision2 = cd2[0];
        } else {
            double cde2 = cd2[oldState];
            double escapeCheck = Math.abs(cde2 - r2);
            double cdc2 = cd2[oldState - 1];
            double coreCheck = Math.abs(cdc2 - r2);
            core = coreCheck < escapeCheck;
            rCollision2 = core ? cdc2 : cde2;
        }
        double ke = bij * bij * reducedMass / (2.0 * rCollision2);
        du[0] = 0;
        int newState = decideBump(atom1, atom2, oldState, core, ke, reducedMass, bij, r2, du, virial, falseTime);
        double lastCollisionVirialr2 = virial[0] / r2;
        Vector dp = Vector.d(v12.getD());
        //dp is the change in momentum due to collision
        dp.Ea1Tv1(lastCollisionVirialr2, r12);
        atom1.getVelocity().PEa1Tv1(rm0, dp);
        atom2.getVelocity().PEa1Tv1(-rm1, dp);
        atom1.getPosition().PEa1Tv1(-falseTime * rm0, dp);
        atom2.getPosition().PEa1Tv1(falseTime * rm1, dp);
        return newState < cd2.length ? newState : -1;
    }

    /**
     * Decides what happens during a collision
     *
     * @param atom1       the first atom in the collision
     * @param atom2       the second atom in the collision
     * @param oldState    the old state of the pair
     * @param core        true if the pair is approach the collision distance from the outside
     * @param ke          kinetic energy of the pair in collision direction
     * @param reducedMass reduced mass for the pair
     * @param bij         rij dot vij
     * @param r2          squared distance at the collision
     * @param du          (out parameter) energy change due to the collision
     * @param virial      (out parameter) virial; used to compute momentum change and pressure
     * @return the new state of the pair
     */
    protected int decideBump(IAtomKinetic atom1, IAtomKinetic atom2, int oldState, boolean core, double ke, double reducedMass, double bij, double r2, double[] du, double[] virial, double falseTime) {
        int newState = oldState + (core ? -1 : +1);
        double uJump = getEnergyForState(atom1, atom2, newState) - getEnergyForState(atom1, atom2, oldState);
        if (ke < uJump) {
            // not enough ke; bounce off core
            virial[0] = 2.0 * reducedMass * bij;
            newState = oldState;
        } else {
            // capture or escape
            virial[0] = reducedMass * (bij + (core ? +1 : -1) * Math.sqrt(bij * bij - 2.0 * r2 * uJump / reducedMass));
            du[0] = uJump;
        }
        return newState;
    }

    public double getRange() {
        return Math.sqrt(collisionDistances2[collisionDistances2.length - 1]);
    }

    @Override
    public double u(double r2) {
        double[] cd2 = collisionDistances2;
        int s = cd2.length;
        for (int i = 0; i < cd2.length; i++) {
            if (cd2[i] > r2) {
                s = i;
                break;
            }
        }
        if (fixOverlap && s == 0) s = 1;
        return getEnergyForState(s);
    }

    @Override
    public double getEnergyForState(int state) {
        return state < 0 || state >= energies.length ? 0 : energies[state];
    }

    protected double[] getEnergies(IAtom atom1, IAtom atom2) {
        return energies;
    }

    public double getEnergyForState(IAtom atom1, IAtom atom2, int state) {
        double[] u = getEnergies(atom1, atom2);
        return state < 0 || state >= u.length ? 0 : u[state];
    }

    @Override
    public double du(double r2) {
        return 0;
    }

    @Override
    public double d2u(double r2) {
        return 0;
    }

}
