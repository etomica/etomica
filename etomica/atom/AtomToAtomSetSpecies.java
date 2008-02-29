package etomica.atom;

import java.io.Serializable;

import etomica.api.IBox;
import etomica.species.ISpecies;

public class AtomToAtomSetSpecies implements AtomToAtomSet, AtomToIndex, Serializable {

    private static final long serialVersionUID = 1L;

    public AtomToAtomSetSpecies(ISpecies species) {
        this.species = species;
    }
    
    public AtomSet getAtomSet(IAtom atom) {
        return moleculeList;
    }
    
    public int getIndex(IAtom atom) {
        return atom.getIndex();
    }
    
    public void setBox(IBox box) {
        moleculeList = box.getMoleculeList(species);
    }

    protected AtomSet moleculeList;
    protected final ISpecies species;
}
