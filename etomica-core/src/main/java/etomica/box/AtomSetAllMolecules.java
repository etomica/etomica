/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;

import java.util.AbstractList;

/**
 * Creates a facade that makes a set of molecule lists look like a single list. This class is
 * configured by calling setMoleculeLists(). This should be called after construction, every time
 * one of the molecule lists in the set is changed by adding or removing a molecule, and when a
 * list is added or removed from the set.
 */
public class AtomSetAllMolecules extends AbstractList<IMolecule> implements IMoleculeList {

    private IMoleculeList[] moleculeLists;
    private int[] moleculeTotals;

    /**
     * Constructs an empty list. Subsequent call to setMoleculeLists() is needed to configure this list.
     */
    public AtomSetAllMolecules() {
        moleculeTotals = new int[1];
    }

    /**
     * @param i specification of the desired molecule.
     * @return a molecule as ordered by the species and then the molecules within the species.
     * @throws IndexOutOfBoundsException if i >= getMoleculeCount() or i < 0
     */
    public IMolecule get(int i) {
        if (i >= size() || i < 0)
            throw new IndexOutOfBoundsException("Index: " + i +
                    ", Number of molecules: " + size());
        int nSpecies = moleculeLists.length;
        if (moleculeTotals[0] > i) {
            return moleculeLists[0].get(i);
        }
        for (int iSpecies = 1; iSpecies < nSpecies; iSpecies++) {
            if (moleculeTotals[iSpecies] > i) {
                return moleculeLists[iSpecies].get(i - moleculeTotals[iSpecies - 1]);
            }
        }
        throw new IllegalStateException("how can this be?!?!?!");
    }

    /**
     * @return total number of molecules in the list.
     */
    public int size() {
        return moleculeTotals[moleculeTotals.length - 1];
    }

    /**
     * Configures this list based on the given array of lists. The ordering of this list is obtained by concatenating
     * the given lists.
     *
     * @param newMoleculeLists lists of molecules that form this list.
     */
    public void setMoleculeLists(IMoleculeList[] newMoleculeLists) {
        moleculeLists = newMoleculeLists;
        if (moleculeTotals.length - 1 != moleculeLists.length) {
            moleculeTotals = new int[moleculeLists.length + 1];
        }
        if (moleculeLists.length == 0) {
            return;
        }
        moleculeTotals[0] = moleculeLists[0].size();
        for (int i = 1; i < moleculeTotals.length - 1; i++) {
            moleculeTotals[i] = moleculeTotals[i - 1] + moleculeLists[i].size();
        }
        moleculeTotals[moleculeTotals.length - 1] = moleculeTotals[moleculeTotals.length - 2];
    }
}
