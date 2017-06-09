/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;

public class AtomSetAllMolecules implements IMoleculeList {

    private static final long serialVersionUID = 1L;
    protected IMoleculeList[] moleculeLists;
    protected int[] moleculeTotals;

    public AtomSetAllMolecules() {
        moleculeTotals = new int[1];
    }

    public IMolecule getMolecule(int i) {
        if (i >= getMoleculeCount() || i < 0)
            throw new IndexOutOfBoundsException("Index: " + i +
                    ", Number of molecules: " + getMoleculeCount());
        int nSpecies = moleculeLists.length;
        if (moleculeTotals[0] > i) {
            return moleculeLists[0].getMolecule(i);
        }
        for (int iSpecies = 1; iSpecies < nSpecies; iSpecies++) {
            if (moleculeTotals[iSpecies] > i) {
                return moleculeLists[iSpecies].getMolecule(i - moleculeTotals[iSpecies - 1]);
            }
        }
        throw new IllegalStateException("how can this be?!?!?!");
    }

    public int getMoleculeCount() {
        return moleculeTotals[moleculeTotals.length - 1];
    }

    public void setMoleculeLists(IMoleculeList[] newMoleculeLists) {
        moleculeLists = newMoleculeLists;
        if (moleculeTotals.length - 1 != moleculeLists.length) {
            moleculeTotals = new int[moleculeLists.length + 1];
        }
        if (moleculeLists.length == 0) {
            return;
        }
        moleculeTotals[0] = moleculeLists[0].getMoleculeCount();
        for (int i = 1; i < moleculeTotals.length - 1; i++) {
            moleculeTotals[i] = moleculeTotals[i - 1] + moleculeLists[i].getMoleculeCount();
        }
        moleculeTotals[moleculeTotals.length - 1] = moleculeTotals[moleculeTotals.length - 2];
    }
}
