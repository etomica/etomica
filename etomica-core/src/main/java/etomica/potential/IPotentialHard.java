/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

import etomica.atom.IAtomKinetic;
import etomica.space.Vector;

/**
 * Interface for a hard pair potential.
 */
public interface IPotentialHard extends IPotentialAtomic {

    /**
     * Returns the state of the pair of atoms (atom1 and atom2) at distance r12
     */
    int getState(IAtomKinetic atom1, IAtomKinetic atom2, Vector r12);

    /**
     * Returns the pair energy for the given state
     */
    double getEnergyForState(int state);

    /**
     * @param atom1 the first atom
     * @param atom2 the second atom
     * @param r12   the position vector between the atoms (r2 - r1)
     * @param v12   the velocity difference between the atoms (v2 - v1)
     * @param state the current state of the pair (overlapped, in well, etc).
     *              the numbering scheme of states is an implementation
     *              detail of the potential class, returned by getState or bump.
     * @return the next collision time for a pair of atoms.
     */
    double collisionTime(IAtomKinetic atom1, IAtomKinetic atom2, Vector r12, Vector v12, int state);

    /**
     * Handles the collision between a pair of atoms.  The atom positions and velocities
     * are updated and the new state is returned.
     *
     * @param atom1     the first atom
     * @param atom2     the second atom
     * @param oldState  the state of the pair before the collision
     * @param r12       the vector between the pair at the collision (r2 - r1)
     * @param v12       the velocity difference between the atoms (v2 - v1)
     * @param falseTime the time into the future when the particles collide
     *                  the particle positions will be updated such that
     *                  r(falsetime) = r(0) + vnew * falsetime
     * @param virial    out parameter for the collision virial
     * @param du        the change in energy from this collision
     * @return the new state of the pair
     */
    int bump(IAtomKinetic atom1, IAtomKinetic atom2, int oldState, Vector r12, Vector v12, double falseTime, double[] virial, double[] du);

    default int nBody() {
        return 2;
    }
}
