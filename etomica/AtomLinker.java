package etomica;

/**
 * Class for constructing linked lists of Atoms.
 * Each Linker points to one atom and another Linker, the next one in the list.
 * Although each atom has built-in ability to link to one next and one previous atom, these
 * Linkers are needed to construct other lists of atoms, particularly for neighbor lists.
 *
 * @author David Kofke
 */
public final class AtomLinker implements java.io.Serializable {
    private final Atom atom;
    private AtomLinker next = null;
    //Constructors
    public AtomLinker(Atom a) {atom = a;}
    public AtomLinker(Atom a, AtomLinker l) {atom = a; next = l;}
    //Access methods
    public Atom atom() {return atom;}
    public AtomLinker next() {return next;}
    public void setNext(AtomLinker l) {next = l;}
}
