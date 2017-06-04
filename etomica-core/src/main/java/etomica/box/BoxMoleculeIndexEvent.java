/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.api.IMolecule;


/**
 * Box event that indicates that a molecule's index has changed.
 */
public class BoxMoleculeIndexEvent extends BoxMoleculeEvent {

    protected final int index;

    public BoxMoleculeIndexEvent(Box box, IMolecule mole, int index) {
        super(box, mole);
        this.index = index;
    }

    /**
     * @return the new index of the molecule
     */
    public int getIndex() {
        return index;
    }
}
