package etomica.box;

import etomica.atom.AtomSet;
import etomica.atom.IAtom;

public class AtomSetAllMoleculesFixed implements AtomSet {

    public AtomSetAllMoleculesFixed() {
    }

    public IAtom getAtom(int i) {
        return allMolecules.getAtom(i);
    }

    public int getAtomCount() {
        return allMolecules.getAtomCount();
    }

    public void setMoleculeLists(AtomSet[] newMoleculeLists) {
        if (newMoleculeLists.length > 0) {
            allMolecules = newMoleculeLists[0];
        }
    }
    
    protected AtomSet allMolecules;
}
