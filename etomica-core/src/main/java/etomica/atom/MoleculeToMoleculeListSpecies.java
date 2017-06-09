/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import java.io.Serializable;

import etomica.box.Box;
import etomica.api.ISpecies;

public class MoleculeToMoleculeListSpecies implements MoleculeToMoleculeList, MoleculeToIndex, Serializable {

    private static final long serialVersionUID = 1L;

    public MoleculeToMoleculeListSpecies(ISpecies species) {
        this.species = species;
    }
    
    public IMoleculeList getMoleculeList(IMolecule molecule) {
        return moleculeList;
    }
    
    public int getIndex(IMolecule atom) {
        return atom.getIndex();
    }
    
    public void setBox(Box box) {
        moleculeList = box.getMoleculeList(species);
    }

    protected IMoleculeList moleculeList;
    protected final ISpecies species;
}
