/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;

import etomica.atom.IAtomList;
import etomica.space.Vector;

/**
 * Returns the position of the first child leaf atom.  Recurses to find
 * the first child leaf atom.
 */

public class MoleculePositionFirstAtom implements IMoleculePositionDefinition, java.io.Serializable {

    public Vector position(IMolecule molecule) {
        IAtomList childList = molecule.getChildList();
        if (childList.getAtomCount() == 0) {
            return null;
        }
        return childList.getAtom(0).getPosition();
    }
    

    private static final long serialVersionUID = 1L;
}
