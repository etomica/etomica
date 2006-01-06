package etomica.atom;

public class AtomToParentChildList implements AtomToArrayList {

    public AtomArrayList getArrayList(Atom atom) {
        return atom.node.parentNode().childList;
    }

}
