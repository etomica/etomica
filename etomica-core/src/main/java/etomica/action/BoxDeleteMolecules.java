/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import etomica.atom.AtomFilter;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;

/**
 * Deletes molecules from a box as determined by an AtomFilter. Atoms deleted
 * are those for which the filter's accept method returns false.
 * 
 * @author David Kofke
 *  
 */
public class BoxDeleteMolecules extends BoxActionAdapter {

    /**
     * @param filter
     *            determines the atoms that will be deleted by the action; those
     *            for which filter.accept returns false are deleted
     */
    public BoxDeleteMolecules(AtomFilter filter) {
        this.filter = filter;
    }

    /**
     * Performs the action of deleting accept == false molecules, considering
     * all molecules in the box last given to setBox. If no box was given,
     * no action is performed and method returns quietly.
     */
    public void actionPerformed() {
        IMoleculeList molecules = box.getMoleculeList();
        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            IMolecule molecule = molecules.getMolecule(i);
            if (!filter.accept(molecule)) {
                box.removeMolecule(molecule);
            }
        }
    }

    private static final long serialVersionUID = 1L;
    private final AtomFilter filter;
} 
