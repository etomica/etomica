package etomica.atom;

public class AtomToParentChildList implements AtomToArrayList, java.io.Serializable {

    public AtomArrayList getArrayList(Atom atom) {
        return atom.node.parentNode().childList;
    }

}
