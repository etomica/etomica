package etomica;

/**
 * Class for constructing linked lists of Atoms.
 *
 * @author David Kofke
 */
public class AtomLinker implements java.io.Serializable {
    public final Atom atom;
    public AtomLinker next = null, previous = null;
    
    /**
     * Constructor throws exception if given atom is null.  Only
     * AtomLink.Tab instances can have a null atom field.
     */
    protected AtomLinker(Atom a) {
        if(a == null && !(this instanceof Tab))
            throw new IllegalArgumentException("Error: cannot create AtomLinker with null atom");
        atom = a; 
    }
    
    //will add a reservoir of linkers, so this is used as constructor for now
    public static AtomLinker makeLinker(Atom a) {
        return new AtomLinker(a);
    }
    
    public static class Tab extends AtomLinker {
        public Tab nextTab, previousTab;
        public Tab() {
            super(null);
        }
    }
        
}//end of AtomLinker
