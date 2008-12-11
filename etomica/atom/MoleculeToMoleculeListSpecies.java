package etomica.atom;

import java.io.Serializable;

import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
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
    
    public void setBox(IBox box) {
        moleculeList = box.getMoleculeList(species);
    }

    protected IMoleculeList moleculeList;
    protected final ISpecies species;
}
