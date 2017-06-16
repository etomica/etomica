/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import java.io.Serializable;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;

/**
 * Static iterator that returns no atoms.
 * @author kofke
 */
public final class AtomIteratorNull implements AtomIterator, Serializable {

    // prevent instantiation.  Consumers should use the INSTANCE field.
    private AtomIteratorNull() {}
    
    public IAtomList next() {return null;}

    public IAtom nextAtom() {return null;}

    public void reset() {}

    public int size() {return 0;}

    public void unset() {}

    public int nBody() {return 1;}
    
    private static final long serialVersionUID = 1L;
    public static final AtomIteratorNull INSTANCE = new AtomIteratorNull();
}
