/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;

import java.util.function.Predicate;

/**
 * Deletes molecules from a box as determined by an AtomFilter. Atoms deleted
 * are those for which the filter's test method returns false.
 * 
 * @author David Kofke
 *  
 */
public class BoxDeleteMolecules extends BoxActionAdapter {

    /**
     * @param filter
     *            determines the atoms that will be deleted by the action; those
     *            for which filter.test returns false are deleted
     */
    public BoxDeleteMolecules(Predicate<IMolecule> filter) {
        this.filter = filter;
    }

    /**
     * Performs the action of deleting test == false molecules, considering
     * all molecules in the box last given to setBox. If no box was given,
     * no action is performed and method returns quietly.
     */
    public void actionPerformed() {
        IMoleculeList molecules = box.getMoleculeList();
        for (int i = 0; i<molecules.size(); i++) {
            IMolecule molecule = molecules.get(i);
            if (!filter.test(molecule)) {
                box.removeMolecule(molecule);
            }
        }
    }

    private final Predicate<IMolecule> filter;
} 
