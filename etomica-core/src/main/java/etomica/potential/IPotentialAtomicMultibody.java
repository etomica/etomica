/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

/**
 * Interface for an atomic non-additive potential.  Although the potential
 * varies with the angles between the atoms, the potential is spherically
 * symmetric in that the atoms do not have orientation and the potential can be
 * computed from the iteratomic distances.
 *
 * @author Andrew Schultz
 */
public interface IPotentialAtomicMultibody extends IPotentialAtomic {
    
    /**
     * Returns the energy for a configuration of atoms having the given pair
     * separations (squared distances).  The pairs are listed in the order
     * 0-1, 0-2... 0-n, 1-2...1-n, 2-3...
     */
    public double energy(double[] r2);
}
