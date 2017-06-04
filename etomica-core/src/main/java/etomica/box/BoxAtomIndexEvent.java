/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.atom.IAtom;
import etomica.api.IBoxAtomIndexEvent;


/**
 * Event that conveys that the maximum global index in a Box has changed.
 */
public class BoxAtomIndexEvent extends BoxAtomEvent implements IBoxAtomIndexEvent {

    public BoxAtomIndexEvent(Box box, IAtom atom, int _index) {
        super(box, atom);
        this.index = _index;
    }

    public int getIndex() {
        return index;
    }
    
    protected int index = -1;
    private static final long serialVersionUID = 1L;
}
