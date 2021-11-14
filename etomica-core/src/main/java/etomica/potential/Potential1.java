/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.space.Space;

/**
 * Potential acting on a single atom or atom group.
 *
 * @author David Kofke
 */
public abstract class Potential1 extends Potential {

    public Potential1(Space space) {
        super(space);
    }

    /**
     * Returns zero.
     */
    public double getRange() {
        return 0.0;
    }
    
    /**
     * Marker interface indicating that a one-body potential is an intramolecular
     * potential, and not, e.g., a potential of interaction with an external field.
     * This is useful when computing energy changes for molecule translations and
     * rotations, for which intramolecular contributions can be ignored.
     */
    public interface Intramolecular {}

}
