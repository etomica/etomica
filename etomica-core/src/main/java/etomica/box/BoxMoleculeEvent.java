/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.api.IMolecule;

/**
 * A box event that is somehow related to a molecule.  The molecule might have
 * been added or removed, or may have a new index.  Details may be determined
 * from the other interfaces implemented by the event object or obtained from
 * calling methods from those interfaces.
 */
public class BoxMoleculeEvent extends BoxEvent {

    private static final long serialVersionUID = 1L;
    protected IMolecule molecule = null;


    public BoxMoleculeEvent(Box box, IMolecule mole) {
        super(box);
        this.molecule = mole;
    }

    /**
     * @return the molecule that is related to this event.
     */
    public IMolecule getMolecule() {
        return molecule;
    }
}
