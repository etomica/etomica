package etomica.atom;

import etomica.api.IAtomLeaf;
import etomica.api.IAtomList;

public class AtomToParentChildList implements AtomToAtomLeafList, java.io.Serializable {

    public IAtomList getAtomList(IAtomLeaf atom) {
        return atom.getParentGroup().getChildList();
    }

    private static final long serialVersionUID = 1L;
}
