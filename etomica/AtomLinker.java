package etomica;

/**
 * Class for constructing linked lists of Atoms.
 * Each Linker points to one atom and another Linker, the next one in the list.
 * Although each atom has built-in ability to link to one next and one previous atom, these
 * Linkers are needed to construct other lists of atoms, particularly for neighbor lists.
 *
 * @author David Kofke
 */
public class AtomLinker implements java.io.Serializable {
    public final Atom atom;
    public AtomLinker next = null, previous = null;
    
    //Constructor
    protected AtomLinker(Atom a) {
        atom = a; 
    }
    
    //will add a reservoir of linkers, so this is used as constructor for now
    public static AtomLinker makeLinker(Atom a) {
        return new AtomLinker(a);
    }
    
    public static class Index extends AtomLinker {
        public Index() {
            super(null);
        }
    }
        
}//end of AtomLinker
