/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.AtomSourceRandomLeaf;

public class AtomSourceRandomMonomer extends AtomSourceRandomLeaf {
    
    public void setAssociationManager(AssociationManager associationManager){
    	this.associationManager = associationManager;
    }

	public IAtom getAtom() {
    	IAtomList atoms = associationManager.getAssociatedAtoms();
    	if (atoms.getAtomCount() == list.getAtomCount()) {//all the atoms are dimer
    		return null;
    	}
    	IAtom atom;
    	do {
    		atom = super.getAtom();
    	}
		while (associationManager.getAssociatedAtoms(atom).getAtomCount()> 0);
        return atom;
	}

	protected AssociationManager associationManager;
}
