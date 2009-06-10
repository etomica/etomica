package etomica.association;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IRandom;
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
