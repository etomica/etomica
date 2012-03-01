package etomica.association;

import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
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
