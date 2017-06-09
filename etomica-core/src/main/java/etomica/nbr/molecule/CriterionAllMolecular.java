/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.nbr.molecule;

import etomica.box.Box;
import etomica.atom.IMolecule;
import etomica.atom.IMoleculeList;

/**
 * Specifies that all atoms pairs are to be considered neighbors.  Should
 * not be used for species in which atoms are being added/removed by integrator.
 */
public class CriterionAllMolecular implements NeighborCriterionMolecular, java.io.Serializable {

    /**
     * 
     */
    private static final long serialVersionUID = 1L;

    /**
     * Always returns false, indicating that neighbor list never needs updating.
     * This is appropriate if atoms are never added to or removed from box,
     * because all atoms are always on neighbor list.
     */
    public boolean needUpdate(IMolecule molecule) {
        return false;
    }

    /**
     * Performs no action.
     */
    public void setBox(Box box) {
    }

    /**
     * Always returns false, indicating that neighbor list never needs updating.
     * This is appropriate if atoms are never added to or removed from box,
     * because all atoms are always on neighbor list.
     */
    public boolean unsafe() {
        return false;
    }

    /**
     * Performs no action.
     */
    public void reset(IMolecule molecule) {
    }

    /**
     * Always returns true, indicating that all atoms pairs are neighbors.
     */
    public boolean accept(IMoleculeList pair) {
        return true;
    }
    
}
