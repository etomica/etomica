/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.molecule.IMoleculeList;

public interface IPotentialMolecular {

    /**
     * Returns the range over which the potential applies.  IAtoms with a
     * greater separation do not interact.
     */
    public double getRange();

    /**
     * Returns the interaction energy between the given atoms.  There might be
     * 0, 1, 2 or more molecules in the IMoleculesList.
     */
    public double energy(IMoleculeList molecules);
}
