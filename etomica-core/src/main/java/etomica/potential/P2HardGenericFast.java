/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

import etomica.atom.IAtomKinetic;
import etomica.space.Vector;

/**
 * Potential class for any hard spherical potential.  This behaves like
 * P2HardGeneric, except that the state is recomputed on the fly internally.
 * This is faster than having the integrator actually store the state, so this
 * a simulation with square-well will run faster with this version of the
 * potential.
 */
public class P2HardGenericFast extends P2HardGeneric {

    public P2HardGenericFast(double[] collisionDistances, double[] energies) {
        super(collisionDistances, energies);
    }

    @Override
    public int getState(IAtomKinetic atom1, IAtomKinetic atom2, Vector r12) {
        // we tell the integrator that our state is -1 (not interacting)
        return -1;
    }

    private int getActualState(Vector r12) {
        return super.getState(null, null, r12);
    }

    @Override
    public double collisionTime(IAtomKinetic atom1, IAtomKinetic atom2, Vector r12, Vector v12, int collisionState) {
        // we need to recompute the state.  fast-forward the distance a little bit since
        // we may be asked to compute a collision time immediately after a simulation
        Vector r12new = Vector.d(r12.getD());
        r12new.E(r12);
        r12new.PEa1Tv1(1e-9, v12);
        collisionState = getActualState(r12new);
        return super.collisionTime(atom1, atom2, r12, v12, collisionState);
    }

    @Override
    public int bump(IAtomKinetic atom1, IAtomKinetic atom2, int oldState, Vector r12, Vector v12, double falseTime, double[] virial, double[] du) {
        // we need to recompute the state.  rewind the distance a little bit,
        Vector r12old = Vector.d(r12.getD());
        r12old.E(r12);
        r12old.PEa1Tv1(-falseTime, v12);
        oldState = getActualState(r12old);
        super.bump(atom1, atom2, oldState, r12, v12, falseTime, virial, du);
        return -1;
    }
}
