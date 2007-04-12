package etomica.atom;

public class AtomToParentChildList implements AtomToArrayList, java.io.Serializable {

    public AtomArrayList getArrayList(IAtom atom) {
        return atom.getParentGroup().getChildList();
    }

    private static final long serialVersionUID = 1L;
}
