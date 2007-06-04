package etomica.atom;

public class AtomToParentChildList implements AtomToAtomSet, java.io.Serializable {

    public AtomSet getAtomSet(IAtom atom) {
        return atom.getParentGroup().getChildList();
    }

    private static final long serialVersionUID = 1L;
}
