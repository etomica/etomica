/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomList;
import etomica.api.IPotentialAtomic;
import etomica.space.Tensor;

/**
 * Interface for hard potentials, having impulsive forces.
 */
public interface PotentialHard extends IPotentialAtomic {

    /**
     * Value of the virial from the most recent collision.
     * 
     * @return double virial value
     */
    public double lastCollisionVirial();

    /**
     * Value of the virial from the most recent collision, decomposed into it
     * tensoral elements.
     * 
     * @return Tensor
     */
    public Tensor lastCollisionVirialTensor();

    /**
     * Implements the collision dynamics. The given atom(s) is assumed to be at
     * the point of collision. This method is called to change their momentum
     * according to the action of the collision. Extensions can be defined to
     * instead implement other, perhaps unphysical changes.
     */
    public void bump(IAtomList atom, double falseTime);

    /**
     * Computes the time of collision of the given atom(s) with the hard
     * potential, assuming no intervening collisions. Usually assumes
     * free-flight between collisions.
     */
    public double collisionTime(IAtomList atom, double falseTime);

    /**
     * returns change in potential energy due to the last collision
     */
    public double energyChange();
}
