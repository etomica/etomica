/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.api.IBox;
import etomica.api.IBoxMoleculeEvent;
import etomica.api.IMolecule;


/**
 * Event that conveys that an Atom has been added to a Box.
 */
public class BoxMoleculeAddedEvent extends BoxMoleculeEvent implements IBoxMoleculeEvent {

    public BoxMoleculeAddedEvent(IBox box, IMolecule mole) {
        super(box, mole);
    }

    private static final long serialVersionUID = 1L;
}
