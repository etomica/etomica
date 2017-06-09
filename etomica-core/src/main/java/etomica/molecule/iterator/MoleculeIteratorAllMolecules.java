/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule.iterator;

import etomica.box.Box;

/**
 * Iterator for all the molecules in a box. Loops over all those atoms that
 * lie just below the species agents in the atom tree hierarchy. To iterate the
 * molecules of just one species, use AtomIteratorMolecule.
 * 
 * @author David Kofke
 * @since 02.02.16
 */

public class MoleculeIteratorAllMolecules extends MoleculeIteratorArrayListSimple 
            implements MoleculeIteratorBoxDependent {

    public MoleculeIteratorAllMolecules() {
        super();
    }

    /**
     * Returns a new iterator ready to iterate over the molecules of the given
     * box.
     */
    public MoleculeIteratorAllMolecules(Box box) {
        this();
        setBox(box);
    }

    /**
     * Sets the box having the molecules to be returned as iterates.
     */
    public void setBox(Box box) {
        setList(box.getMoleculeList());
    }
    
    private static final long serialVersionUID = 1L;

}
