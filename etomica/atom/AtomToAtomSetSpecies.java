package etomica.atom;

import java.io.Serializable;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.ISpecies;

public class AtomToAtomSetSpecies implements AtomToAtomSet, AtomToIndex, Serializable {

    private static final long serialVersionUID = 1L;

    public AtomToAtomSetSpecies(ISpecies species) {
        this.species = species;
    }
    
    public IAtomList getAtomSet(IAtom atom) {
        return moleculeList;
    }
    
    public int getIndex(IAtom atom) {
        return ((IMolecule)atom).getIndex();
    }
    
    public void setBox(IBox box) {
        moleculeList = box.getMoleculeList(species);
    }

    protected IAtomList moleculeList;
    protected final ISpecies species;
}
