/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.lattice.CellLattice;

/**
 * Interface for class that handles assignment of atoms to cells in a box.
 * This facility is needed by neighbor-listing schemes, though it may find use for
 * other purposes.
 *
 * @author David Kofke and Andrew Schultz
 */
public interface BoxCellManager {

    /**
     * Returns the lattice that defines the cell arrangement.
     */
    CellLattice getLattice();

    /**
     * Assigns cells to all molecules in the box.
     */
    void assignCellAll();
}
