/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.atom.IAtom;

/**
 * A box event that is somehow related to an atom.  The atom might have been
 * added or removed, or may have a new index.  Details may be determined from
 * the other interfaces implemented by the event object or obtained from
 * calling methods from those interfaces.
 */
public class BoxAtomEvent extends BoxEvent {

    protected final IAtom atom;

    public BoxAtomEvent(Box box, IAtom atom) {
        super(box);
        this.atom = atom;
    }

    /**
     * @return the atom that is related to this event.
     */
    public IAtom getAtom() {
        return atom;
    }
}
