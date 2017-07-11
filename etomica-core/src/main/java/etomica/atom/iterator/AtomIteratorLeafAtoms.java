/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.box.Box;

/**
 * Iterator that will loop over all leaf atoms in a box. Can be configured to
 * iterate all leaf atoms, or only those of a particular species.
 */
public class AtomIteratorLeafAtoms extends AtomIteratorArrayListSimple implements 
        AtomIteratorBoxDependent {

    /**
     * Creates iterator with no box specified. Iteration will return no atoms
     * until a call to setBox is performed.
     */
    public AtomIteratorLeafAtoms() {
        super();
    }

    /**
     * Creates iterator conditioned to give all leaf atoms of the specified
     * box. Call to reset() is required before beginning iteration.
     */
    public AtomIteratorLeafAtoms(Box box) {
        this();
        setBox(box);
    }

    /**
     * Configures iterator to form its iterates from the leaf atoms of the given
     * box. If a species was previously (or subsequently) set, iterates will
     * be the leaf atoms of under the species in the specified box.
     * @throws NullPointerException if the Box is null
     */
    public void setBox(Box box) {
        setList(box.getLeafList());
    }

    private static final long serialVersionUID = 1L;
}
