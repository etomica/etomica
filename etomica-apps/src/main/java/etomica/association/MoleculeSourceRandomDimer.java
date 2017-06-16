/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeSource;
import etomica.util.random.IRandom;

public class MoleculeSourceRandomDimer implements MoleculeSource {
	public void setRandomNumberGenerator(IRandom newRandom) {
        random = newRandom;
    }
    
    /**
     * Returns the random number generator used to pick atoms
     */
    public IRandom getRandomNumberGenerator() {
        return random;
    }
    
    public void setAssociationManager(AssociationManagerMolecule associationManager){
    	this.associationManager = associationManager;
    }

	public IMolecule getMolecule() {
    	IMoleculeList molecules = associationManager.getAssociatedMolecules();
    	if (molecules.getMoleculeCount() == 0) {
    		return null;
    	}
        return molecules.getMolecule(random.nextInt(molecules.getMoleculeCount()));
	}

	public void setBox(Box p) {

	}
	protected IRandom random;
	protected AssociationManagerMolecule associationManager;
}
