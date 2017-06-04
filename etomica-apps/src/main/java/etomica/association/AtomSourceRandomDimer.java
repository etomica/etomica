/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.api.IRandom;
import etomica.atom.AtomSource;

public class AtomSourceRandomDimer implements AtomSource {
	public void setRandomNumberGenerator(IRandom newRandom) {
        random = newRandom;
    }
    
    /**
     * Returns the random number generator used to pick atoms
     */
    public IRandom getRandomNumberGenerator() {
        return random;
    }
    
    public void setAssociationManager(AssociationManager associationManager){
    	this.associationManager = associationManager;
    }

	public IAtom getAtom() {
    	IAtomList atoms = associationManager.getAssociatedAtoms();
    	if (atoms.getAtomCount() == 0) {
    		return null;
    	}
        return atoms.getAtom(random.nextInt(atoms.getAtomCount()));
	}

	public void setBox(Box p) {

	}
	protected IRandom random;
	protected AssociationManager associationManager;
}
