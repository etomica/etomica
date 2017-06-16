/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import java.io.Serializable;

public class AtomToAtomSetFixed implements AtomToAtomLeafList, AtomToIndex, Serializable {

    private static final long serialVersionUID = 1L;

    public AtomToAtomSetFixed() {
        atomArrayList = new AtomArrayList();
    }
    
    public void setArrayList(AtomArrayList list) {
        atomArrayList = list;
    }
    
    public IAtomList getAtomList(IAtom atom) {
        return atomArrayList;
    }
    
    public int getIndex(IAtom atom) {
        return atomArrayList.indexOf(atom);
    }

    private AtomArrayList atomArrayList;
}
