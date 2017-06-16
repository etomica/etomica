/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.box.Box;

/**
 * Iterator that returns all pairs that can be formed from all leaf atoms of a
 * given box.
 */
public class ApiLeafAtoms extends ApiIntraArrayList implements
        AtomsetIteratorBoxDependent {

    /**
     * Creates new pair iterator that requires reset() before beginning
     * iteration.
     */
    public ApiLeafAtoms() {
        super();
    }

    /**
     * Conditions iterator to return all leaf-atom pairs from the given box.
     * @throws NullPointerException if the Box is null
     */
    public void setBox(Box box) {
        setList(box.getLeafList());
    }

    private static final long serialVersionUID = 1L;
}
