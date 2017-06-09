/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.atom.IMolecule;
import etomica.atom.IMoleculeList;
import etomica.atom.MoleculeSourceRandomMolecule;

public class MoleculeSourceRandomMonomer extends MoleculeSourceRandomMolecule {
    
    public void setAssociationManager(AssociationManagerMolecule associationManager){
    	this.associationManager = associationManager;
    }

	public IMolecule getMolecule() {
    	IMoleculeList molecules = associationManager.getAssociatedMolecules();
    	if (molecules.getMoleculeCount() == box.getMoleculeList().getMoleculeCount()) {//all the molecules are dimer
    		return null;
    	}
    	IMolecule molecule;
    	do {
    		molecule = super.getMolecule();
    	}
		while (associationManager.getAssociatedMolecules(molecule).getMoleculeCount()> 0);
        return molecule;
	}

	protected AssociationManagerMolecule associationManager;
}
