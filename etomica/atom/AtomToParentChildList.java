package etomica.atom;

import etomica.api.IAtom;
import etomica.api.IAtomLeaf;
import etomica.api.IAtomList;

public class AtomToParentChildList implements AtomToAtomSet, java.io.Serializable {

    public IAtomList getAtomSet(IAtom atom) {
        return ((IAtomLeaf)atom).getParentGroup().getChildList();
    }

    private static final long serialVersionUID = 1L;
}
