/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.api.IAtom;
import etomica.api.IAtomList;

public class AtomToParentChildList implements AtomToAtomLeafList, java.io.Serializable {

    public IAtomList getAtomList(IAtom atom) {
        return atom.getParentGroup().getChildList();
    }

    private static final long serialVersionUID = 1L;
}
