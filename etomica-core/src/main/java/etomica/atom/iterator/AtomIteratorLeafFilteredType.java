/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IAtomType;
import etomica.api.IBox;

public class AtomIteratorLeafFilteredType extends AtomIteratorLeafAtoms {

    public AtomIteratorLeafFilteredType(IBox box, IAtomType type) {
        super(box);
        filteredType = type;
    }
    
    public IAtom nextAtom() {
        IAtom atom = super.nextAtom();
        while (atom != null) {
            if (atom.getType() == filteredType) {
                return atom;
            }
            atom = super.nextAtom();
        }
        return null;
    }

    public IAtomList next() {
        IAtomList atom = super.next();
        while (atom != null) {
            if (atom.getAtom(0).getType() == filteredType) {
                return atom;
            }
            atom = super.next();
        }
        return null;
    }

    public int size() {
        reset();
        int count = 0;
        for (IAtom atom = nextAtom(); atom != null; atom = nextAtom()) {
            count++;
        }
        return count;
    }

    private static final long serialVersionUID = 1L;
    protected final IAtomType filteredType;
}
