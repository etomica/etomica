/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

import etomica.atom.IAtomKinetic;
import etomica.space.Vector;

/**
 * Interface for a hard external field potential.  The additional methods here
 * facilitate hard dynamics.  The state for a given atom is an implementation
 * detail of the potential, except that a negative state indicates that an atom
 * is not interacting (energy = 0).
 */
public interface IPotentialHardField extends IPotentialField {

    /**
     * Returns the state of the atom
     */
    int getState(IAtomKinetic atom);

    /**
     * Returns the pair energy for the given state
     */
    double getEnergyForState(int state);

    /**
     * Computes the collision time for a given atom.
     *
     * @param atom  the atom
     * @param state the current state of the atom (overlapped, in well, etc).
     *              the numbering scheme of states is an implementation
     *              detail of the potential class, returned by getState or bump.
     * @return the next collision time for a pair of atoms.
     */
    double collisionTime(IAtomKinetic atom, Vector r, Vector v, int state);

    /**
     * Handles the collision between a pair of atoms.  The atom positions and velocities
     * are updated and the new state is returned.
     *
     * @param atom      the atom
     * @param oldState  the state of the pair before the collision
     * @param falseTime the time into the future when the atoms collides
     *                  the particle position will be updated such that
     *                  r(falsetime) = r(0) + vnew * falsetime
     * @param du        the change in energy from this collision
     * @return the new state of the atom
     */
    int bump(IAtomKinetic atom, int oldState, Vector r, double falseTime, Vector deltaP, double[] du);

    default int nBody() {
        return 1;
    }
}
