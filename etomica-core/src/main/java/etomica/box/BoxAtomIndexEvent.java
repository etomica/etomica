/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.atom.IAtom;


/**
 * Box event that indicates that an atom's index has changed.
 */
public class BoxAtomIndexEvent extends BoxAtomEvent {

    protected final int index;

    public BoxAtomIndexEvent(Box box, IAtom atom, int index) {
        super(box, atom);
        this.index = index;
    }

    /**
     * @return the index of the atom
     */
    public int getIndex() {
        return index;
    }
}
