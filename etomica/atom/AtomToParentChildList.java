package etomica.atom;

import etomica.api.IAtom;
import etomica.api.IAtomLeaf;
import etomica.api.IAtomSet;

public class AtomToParentChildList implements AtomToAtomSet, java.io.Serializable {

    public IAtomSet getAtomSet(IAtom atom) {
        return ((IAtomLeaf)atom).getParentGroup().getChildList();
    }

    private static final long serialVersionUID = 1L;
}
