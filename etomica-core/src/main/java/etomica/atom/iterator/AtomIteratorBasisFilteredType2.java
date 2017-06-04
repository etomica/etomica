/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.IAtomType;

/**
 * This iterator returns atoms that are one of two types.  This can be used to
 * construct an ApiIntergroupIntraSpecies iterator.
 * 
 * @author Andrew Schultz
 */
public class AtomIteratorBasisFilteredType2 extends AtomIteratorBasis {

    public AtomIteratorBasisFilteredType2(IAtomType type1, IAtomType type2) {
        super();
        filteredType1 = type1;
        filteredType2 = type2;
    }
    
    public IAtom nextAtom() {
        IAtom atom = super.nextAtom();
        while (atom != null) {
            if (atom.getType() == filteredType1 || atom.getType() == filteredType2) {
                return atom;
            }
            atom = super.nextAtom();
        }
        return null;
    }

    public IAtomList next() {
        IAtomList atom = super.next();
        while (atom != null) {
            IAtomType t = atom.getAtom(0).getType();
            if (t == filteredType1 || t == filteredType2) {
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
    protected final IAtomType filteredType1, filteredType2;
}
