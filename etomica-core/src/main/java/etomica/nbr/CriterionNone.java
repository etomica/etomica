/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;

public final class CriterionNone implements NeighborCriterion, java.io.Serializable {
    private static final long serialVersionUID = 1L;

    /**
     * Always returns false, indicating that neighbor list never needs updating.
     * This is appropriate if atoms are never added to or removed from box,
     * because all atoms are always on neighbor list.
     */
    public boolean needUpdate(IAtom atom) {return false;}

    /**
     * Performs no action.
     */
    public void setBox(Box box) {}

    /**
     * Always returns false, indicating that neighbor list never needs updating.
     */
    public boolean unsafe() {return false;}

    /**
     * Performs no action.
     */
    public void reset(IAtom atom) {}

    /**
     * Always returns false, indicating that no atoms pairs are neighbors.
     */
    public boolean accept(IAtomList pair) {return false;}
}
