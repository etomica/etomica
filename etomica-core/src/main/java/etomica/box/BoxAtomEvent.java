/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.atom.IAtom;
import etomica.api.IBoxAtomEvent;

/**
 * Event that conveys some happening with respect to an Atom in a Box.
 */
public class BoxAtomEvent extends BoxEvent implements IBoxAtomEvent {
    
    public BoxAtomEvent(Box box, IAtom atom) {
        super(box);
        this.atom = atom;
    }

    public IAtom getAtom() {
        return atom;
    }
    
    protected IAtom atom = null;
    private static final long serialVersionUID = 1L;
}
