package etomica.box;

import etomica.atom.AtomSet;
import etomica.atom.IAtom;

public class AtomSetAllMolecules implements AtomSet {

    public AtomSetAllMolecules() {
        moleculeTotals = new int[1];
    }

    public IAtom getAtom(int i) {
        if(i >= getAtomCount() || i < 0) 
            throw new IndexOutOfBoundsException("Index: "+i+
                                                ", Number of molecules: "+getAtomCount());
        int nSpecies = moleculeLists.length;
        if (moleculeTotals[0] > i) {
            return moleculeLists[0].getAtom(i);
        }
        for (int iSpecies=1; iSpecies<nSpecies; iSpecies++) {
            if (moleculeTotals[iSpecies] > i) {
                return moleculeLists[iSpecies].getAtom(i-moleculeTotals[iSpecies-1]);
            }
        }
        throw new IllegalStateException("how can this be?!?!?!");
    }

    public int getAtomCount() {
        return moleculeTotals[moleculeTotals.length-1];
    }

    public void setMoleculeLists(AtomSet[] newMoleculeLists) {
        moleculeLists = newMoleculeLists;
        if (moleculeTotals.length - 1 != moleculeLists.length) {
            moleculeTotals = new int[moleculeLists.length+1];
        }
        if (moleculeLists.length == 0) {
            return;
        }
        moleculeTotals[0] = moleculeLists[0].getAtomCount();
        for (int i=1; i<moleculeTotals.length-1; i++) {
            moleculeTotals[i] = moleculeTotals[i-1] + moleculeLists[i].getAtomCount();
        }
        moleculeTotals[moleculeTotals.length-1] = moleculeTotals[moleculeTotals.length-2];
    }
    
    protected AtomSet[] moleculeLists;
    protected int[] moleculeTotals;
}
