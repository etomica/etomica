package etomica;

/**
 * Class for making linked lists of AtomPairs
 */
public final class AtomPairLinker implements java.io.Serializable {
    private AtomPair pair;
    private AtomPairLinker next;
    public AtomPairLinker() {}
    public AtomPairLinker(AtomPair p) {pair = p;}
    public AtomPairLinker(AtomPair p, AtomPairLinker l) {pair = p; next = l;}
    public final AtomPair pair() {return pair;}
    public final AtomPairLinker next() {return next;}
    public final void setNext(AtomPairLinker l) {next = l;}
    public final void setPair(AtomPair p) {pair = p;}
} //end of Linker
