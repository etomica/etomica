/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.space.Vector;

public interface IPotentialField {

    /**
     * Returns the energy between IAtom atom and the field.
     */
    double u(IAtom atom);

    /**
     * Computes the force (and adds it to f) for IAtom atom and returns the
     * energy due to the field.
     */
    double udu(IAtom atom, Vector f);

    default double uduTorque(IAtom atom, Vector f, Vector t) {
        return udu(atom, f);
    }

    /**
     * Returns the state of the atom
     */
    default int getState(IAtomKinetic atom) {
        throw new RuntimeException("Atoms have a state only for hard potentials");
    }

    /**
     * Returns the pair energy for the given state
     */
    default double getEnergyForState(int state) {
        throw new RuntimeException("Atoms have a state only for hard potentials");
    }

    /**
     * Computes the collision time for a given atom.
     *
     * @param atom      the atom
     * @param state     the current state of the atom (overlapped, in well, etc).
     *                  the numbering scheme of states is an implementation
     *                  detail of the potential class, returned by getState or bump.
     * @param falseTime
     * @return the next collision time for a pair of atoms.
     */
    default double collisionTime(IAtomKinetic atom, Vector r, Vector v, int state, double falseTime) {
        throw new RuntimeException("Atoms have a collision time only for hard potentials");
    }

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
    default int bump(IAtomKinetic atom, int oldState, Vector r, double falseTime, Vector deltaP, double[] du) {
        throw new RuntimeException("Atoms can bump only for hard potentials");
    }

}
