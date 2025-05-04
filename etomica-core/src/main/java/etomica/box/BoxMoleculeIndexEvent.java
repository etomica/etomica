/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.molecule.IMolecule;


/**
 * Box event that indicates that a molecule's index has changed.
 */
public class BoxMoleculeIndexEvent extends BoxMoleculeEvent {

    protected final int oldIndex;

    public BoxMoleculeIndexEvent(Box box, IMolecule mole, int oldIndex) {
        super(box, mole);
        this.oldIndex = oldIndex;
    }

    /**
     * @return the new index of the molecule
     */
    public int getOldIndex() {
        return oldIndex;
    }
}
