/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential.compute;

/**
 * Interface for neighbor manager for hard MD.
 */
public interface NeighborManagerHard extends NeighborManager {

    /**
     * Updates the state of pair of atoms i and j.
     */
    void setPairState(int i, int j, int state);
}
