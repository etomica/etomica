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
    public final Atom atom;
    public AtomLinker next = null, previous = null;
    //Constructor
    private AtomLinker(Atom a, AtomLinker next, AtomLinker previous) {
        atom = a; 
        this.next = next;
        this.previous = previous;
        if(next != null) next.previous = this;
        if(previous != null) previous.next = this;
    }
    
    //will add a reservoir of linkers, so this is used as constructor for now
    public static AtomLinker makeLinker(Atom a, AtomLinker next, AtomLinker previous) {
        return new AtomLinker(a, next, previous);
    }
    
    public void delete() {
        if(previous != null) previous.next = next;
        if(next != null) next.previous = previous;
    }
}
