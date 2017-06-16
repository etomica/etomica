/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

/**
 * This provides an interface for a class that acts like a hashmap for atomic
 * diameters.
 *
 * @author Andrew Schultz
 */
public interface DiameterHash {

    /**
     * Returns an appropriate diameter for the given atom.  If no diameter is
     * known for the atom, -1 is returned.
     */
    public double getDiameter(IAtom atom);

}
