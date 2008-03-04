package etomica.box;

import etomica.api.IAtom;
import etomica.api.IAtomSet;
import etomica.atom.AtomSet;

public class AtomSetAllMoleculesFixed implements AtomSet {

    public AtomSetAllMoleculesFixed() {
    }

    public IAtom getAtom(int i) {
        return allMolecules.getAtom(i);
    }

    public int getAtomCount() {
        return allMolecules.getAtomCount();
    }

    public void setMoleculeLists(IAtomSet[] newMoleculeLists) {
        if (newMoleculeLists.length > 0) {
            allMolecules = newMoleculeLists[0];
        }
    }
    
    protected IAtomSet allMolecules;
}
